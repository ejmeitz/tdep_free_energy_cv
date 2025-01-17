import numpy as np
import logging
import os, shutil
from os.path import join

from src import (
    ExtractForceConstants,
    LammpsSimulator,
    get_n_atoms_from_dump,
    remove_dump_headers,
    write_tdep_meta
)

from .configs import InterpolateIFCParams

def run_lammps(p, T, base_infile_path, sim_root_dir):
    N_steps = p.lds.n_configs * p.lds.data_interval
    var_dict = {"T" : T, "N_steps" : N_steps, "structure_path" : p.lds.structure_path}
    required_vars = ["T", "N_steps", "structure_path"]
    ls = LammpsSimulator(base_infile_path, sim_root_dir, var_dict, required_vars)

    ls.mpirun(sim_root_dir, p.lds.n_cores)
    N_atoms = get_n_atoms_from_dump(join(sim_root_dir, "dump.positions"))

    remove_dump_headers(sim_root_dir) # makes infile.positions, infile.stat, infile.forces
    write_tdep_meta(sim_root_dir, N_atoms, p.lds.n_configs, p.lds.time_step_fs, T)

    # Parse thermodata.txt and get actual average temp
    temps = np.loadtxt(join(sim_root_dir, "thermo_data.txt"))[: , p.lds.thermo_data_temp_idx]
    return np.mean(temps)

def write_temps_file(out_dir, temps_to_calculate, mean_sim_temps):
    with open(join(out_dir, "ACTUAL_SIM_TEMPS.txt"), "w") as f:
        f.write("# Requested Temperatures:\n")
        for T in temps_to_calculate:
            f.write(f"{T} ")
        f.write("\n")
        f.write("# Mean Simulation Temperatures:\n")
        for T in mean_sim_temps:
            f.write(f"{T} ")
        f.write("\n")

def run_interpolate_irred(p : InterpolateIFCParams):

    # Check that interp temps are not actually extrapolating
    min_T_calc = np.amin(p.temps_to_calculate)
    max_T_calc = np.amax(p.temps_to_calculate)
    min_T_interp = np.amin(p.temps_to_interpolate)
    max_T_interp = np.amax(p.temps_to_interpolate)

    if max_T_interp > max_T_calc or min_T_interp < min_T_calc:
        raise ValueError("Temperatures passed require extrapolation, exiting.")
    
    TEMPS_CALC = np.sort(p.temps_to_calculate)
    TEMPS_INTERP = np.sort(p.temps_to_interpolate)

    
    # Pre-calc stuff needed for force calculators
    if p.force_calc == "lammps":
        if p.lds is None:
            raise ValueError("Force calculator set to lammps, must pass LammpsDynamicsSettings")
        file_path = os.path.dirname(os.path.realpath(__file__))
        base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.lds.script_name)
        base_infile_path = os.path.abspath(base_infile_path)
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")

    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4)

    sim_dir = join(p.basepath, "simulation_data")
    os.mkdir(sim_dir)

    mean_sim_temps = []
    irred2_ifcs = []
    irred3_ifcs = []
    irred4_ifcs = []
    U0s = []
    for i, T in enumerate(TEMPS_CALC):
        sim_root_dir = join(sim_dir, f"T{T}")
        os.mkdir(sim_root_dir)

        # Move files needed for extract IFCs
        shutil.copyfile(join(p.basepath, "infile.ucposcar"), 
                        join(sim_root_dir, "infile.ucposcar"))
        shutil.copyfile(join(p.basepath, "infile.ssposcar"), 
                        join(sim_root_dir, "infile.ssposcar"))

        if p.force_calc == "lammps":
            logging.info(f"Running LAMMPS at {T} K")
            mean_sim_temp = run_lammps(p, T, base_infile_path, sim_root_dir)
            mean_sim_temps.append(mean_sim_temp)
        else:
            raise NotImplementedError()

        # Calculate force constants from lammps simulation
        res = ef.mpirun(p.n_cores, sim_root_dir)
        if res != 0:
            logging.error(f"Possible error when extracting IFCs")

        # Store irred ifcs
        tmp = np.loadtxt(join(sim_root_dir, f"outfile.irrifc_secondorder"))
        irred2_ifcs.append(tmp)

        if p.rc3 is not None:
            tmp = np.loadtxt(join(sim_root_dir, f"outfile.irrifc_thirdorder"))
            irred3_ifcs.append(tmp)
        
        if p.rc4 is not None:
            tmp = np.loadtxt(join(sim_root_dir, f"outfile.irrifc_fourthorder"))
            irred4_ifcs.append(tmp)

        if p.interpolate_U0:
            tmp = np.loadtxt(join(sim_root_dir, "outfile.U0"))
            U0s.append(tmp[-1])

    # Save actual temps to see if our IFCs are accurate
    if p.force_calc == "lammps":       
        write_temps_file(p.basepath, TEMPS_CALC, mean_sim_temps)
        logging.info("Using mean simulation temperatures for interpolation.")

    # Interpolate the irreducible IFCs
    irred2_ifcs = np.array(irred2_ifcs)
    new_irred2_ifcs = np.zeros((irred2_ifcs.shape[1], len(TEMPS_INTERP)))
    for i in range(irred2_ifcs.shape[1]):
        new_irred2_ifcs[i,:] = np.interp(TEMPS_INTERP, mean_sim_temps, irred2_ifcs[:,i])

    if p.rc3 is not None:
        irred3_ifcs = np.array(irred3_ifcs)
        new_irred3_ifcs = np.zeros((irred3_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred2_ifcs.shape[1]):
            new_irred3_ifcs[i,:] = np.interp(TEMPS_INTERP, mean_sim_temps, irred3_ifcs[:,i])

    if p.rc4 is not None:
        irred4_ifcs = np.array(irred4_ifcs)
        new_irred4_ifcs = np.zeros((irred4_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred2_ifcs.shape[1]):
            new_irred4_ifcs[i,:] = np.interp(TEMPS_INTERP, mean_sim_temps, irred4_ifcs[:,i])

    if p.interpolate_U0:
        interpolated_U0s = np.interp(TEMPS_INTERP, mean_sim_temps, U0s)

    # Save interpolated IFCs per tempearture
    irred2_paths = []
    irred3_paths = [] if p.rc3 is not None else None
    irred4_paths = [] if p.rc4 is not None else None
    irred_out_path = join(p.basepath, "INTERPOLATED_IRRED_IFCS")
    os.mkdir(irred_out_path)
    for j, T in enumerate(TEMPS_INTERP):
        temp_str = f"{T}".replace('.', '_')
        p2 = join(irred_out_path, f"infile.irrifc_secondorder_{temp_str}")
        np.savetxt(p2, new_irred2_ifcs[:,j])
        irred2_paths.append(p2)

        if p.rc3 is not None:
            p3 = join(irred_out_path, f"infile.irrifc_thirdorder_{temp_str}")
            np.savetxt(p3, new_irred3_ifcs[:,j])
            irred3_paths.append(p3)

        if p.rc4 is not None:
            p4 = join(irred_out_path, f"infile.irrifc_fourthorder_{temp_str}")
            np.savetxt(p4, new_irred4_ifcs[:,j])
            irred4_paths.append(p4)
    
    if p.interpolate_U0:
        np.savetxt(join(irred_out_path, "interpolated_U0s.txt"), np.column_stack([TEMPS_INTERP,interpolated_U0s]),
                     fmt = "%.10f", header = "# Temperature U0")

    # WRITE EXTRACT IFC COMMAND USED 
    # THIS MUST BE THE SAME WHEN USING 
    # THE IRRED IFC CREATED HERE
    ef_tmp = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    with open(join(p.basepath, "READ_IRRED_IFC_CMD.txt"), "w") as f:
        f.write("""
                # YOU MUST USE THIS COMMAND WHEN USING THE INTERPOLATED IFCS.\n
                # THE PARAMS MUST MATCH THOSE USED TO CREATE THE INTERPOLATION POINTS.\n
                # THE infile.ucposcar AND infile.ssposcar MUST ALSO BE IDENTICAL.\n"""
                )
        f.write(ef_tmp._cmd())


    if p.make_ss_ifcs:
        logging.info("Converting interpolated irreducible IFCs, back to super cell IFCs")
        ss_out_dir = convert_irred_to_ss_ifc(p, TEMPS_INTERP, TEMPS_CALC, sim_dir, irred_out_path)
        shutil.copyfile(join(irred_out_path, "interpolated_U0s.txt"), join(ss_out_dir, "interpolated_U0s.txt"))

    if p.cleanup:
        logging.info("Cleaning up...will delete simulation data.")
        os.system(f"rm -rf {sim_dir}")

    return irred_out_path, ss_out_dir
    
def convert_irred_to_ss_ifc(p : InterpolateIFCParams, TEMPS_INTERP, TEMPS_CALC, sim_dir, irred_out_path):
    
    ss_out_dir = join(p.basepath, "INTERPOLATED_SUPERCELL_IFCS")
    os.mkdir(ss_out_dir)

    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    for T in TEMPS_INTERP:
        temp_str = f"{T}".replace('.', '_')
        sim_dir_T = join(sim_dir, f"T{temp_str}_SS_RECONSTRUCT")
        os.mkdir(sim_dir_T)
        os.chdir(sim_dir_T)

        # Copy interpolated irreducible IFCs to this dir
        shutil.copyfile(join(irred_out_path, f"infile.irrifc_secondorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_secondorder"))
        if p.rc3 is not None:
            shutil.copyfile(join(irred_out_path, f"infile.irrifc_thirdorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_thirdorder"))
        if p.rc4 is not None:
            shutil.copyfile(join(irred_out_path, f"infile.irrifc_fourthorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_fourthorder"))

        # Move required files
        shutil.copyfile(join(p.basepath, "infile.ucposcar"), join(sim_dir_T, "infile.ucposcar"))
        shutil.copyfile(join(p.basepath, "infile.ssposcar"), join(sim_dir_T, "infile.ssposcar"))

        # forces, positions and stat file will not affect
        # IFC calculation but need to be valid. Just
        # use the files from the first temp
        shutil.copyfile(join(sim_dir, f"T{TEMPS_CALC[0]}", "infile.positions"), join(sim_dir_T, "infile.positions"))
        shutil.copyfile(join(sim_dir, f"T{TEMPS_CALC[0]}", "infile.forces"), join(sim_dir_T, "infile.forces"))
        stat_fmt = ['%d', '%d'] + ['%.8f'] * 11
        np.savetxt(join(sim_dir_T, "infile.stat"), np.zeros((p.lds.n_configs, 13)), fmt = stat_fmt) # this isnt even used, but still has to exist 

        # Write a new meta file cause we can
        write_tdep_meta(sim_dir_T, p.lds.n_atoms, p.lds.n_configs, p.lds.time_step_fs, T)

        # Re-run extract force constants using the irreducible IFCs
        # This still needs infile.meta, infile.stat, infile.positions and infile.forces
        # to run although those will not be used for calculating the IFCs
        ef.mpirun(p.n_cores, sim_dir_T)

        # Move force constants to make folder to avoid clean up 
        shutil.copyfile(join(sim_dir_T, "outfile.forceconstant"), join(ss_out_dir, f"infile.forceconstant"))
        
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_thirdorder"), join(ss_out_dir, f"infile.forceconstant_thirdorder"))
            
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_thirdorder"), join(ss_out_dir, f"infile.forceconstant_thirdorder"))

    return ss_out_dir
