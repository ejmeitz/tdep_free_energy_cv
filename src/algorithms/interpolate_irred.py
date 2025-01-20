import numpy as np
from numpy.polynomial.polynomial import Polynomial
import scipy.interpolate as interp
import logging
import os, shutil
from os.path import join

from src import (
    ExtractForceConstants,
    run_ifc_from_MD,
    IFC_MD_Params,
    LammpsSimulator,
    get_n_atoms_from_dump,
    remove_dump_headers,
    write_tdep_meta,
    temp_to_str
)

from .configs import InterpolateIFCParams, Paths


def interpolate(mode, X_INTERP, X_DATA, Y_DATA):

    if mode == "linear":
        return np.interp(X_INTERP, X_DATA, Y_DATA)
    elif mode == "lagrange":
        poly = Polynomial(interp.lagrange(X_DATA, Y_DATA).coef[::-1])
        return poly(X_INTERP)
    else:
        raise ValueError(f"Invalid interpolate mode : {mode}. Expected linear or lagrange.")

def run_interpolate_irred(p : InterpolateIFCParams, paths : Paths):

    # Check that interp temps are not actually extrapolating
    min_T_calc = np.amin(p.temps_to_simulate)
    max_T_calc = np.amax(p.temps_to_simulate)
    min_T_interp = np.amin(p.temps_to_interpolate)
    max_T_interp = np.amax(p.temps_to_interpolate)

    if max_T_interp > max_T_calc or min_T_interp < min_T_calc:
        raise ValueError("Temperatures passed require extrapolation, exiting.")
    
    TEMPS_SIM = np.sort(p.temps_to_simulate)
    TEMPS_INTERP = np.sort(p.temps_to_interpolate)

    # Get force constants
    if p.force_calc == "lammps":
        # We will clean up manually at end, might need this sim data later
        ifc_params = IFC_MD_Params(p.temps_to_simulate, p.n_cores_max, p.rc2, p.rc3, p.rc4, p.lds, False)
        sim_dir, ifc_dir, mean_sim_temps = run_ifc_from_MD(ifc_params, paths)
        logging.info("Using mean simulation temperatures for interpolation.")
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")

    # Parse irreducible force constants
    irred2_ifcs = []; irred3_ifcs = []; irred4_ifcs = []; U0s = []
    N_atoms = None
    for i, T in enumerate(TEMPS_SIM):
        temp_str = temp_to_str(T)
        ifc_dir_T = join(ifc_dir, f"T{temp_str}")

        # Store irred ifcs
        irred2_ifcs.append(np.loadtxt(join(ifc_dir_T, f"outfile.irrifc_secondorder")))

        if p.rc3 is not None:
            irred3_ifcs.append(tmp = np.loadtxt(join(ifc_dir_T, f"outfile.irrifc_thirdorder")))
        
        if p.rc4 is not None:
            irred4_ifcs.append(np.loadtxt(join(ifc_dir_T, f"outfile.irrifc_fourthorder")))

        if p.interpolate_U0:
            tmp = np.loadtxt(join(ifc_dir_T, "outfile.U0"))
            U0s.append(tmp[-1])

    # Interpolate the irreducible IFCs
    irred2_ifcs = np.array(irred2_ifcs)
    new_irred2_ifcs = np.zeros((irred2_ifcs.shape[1], len(TEMPS_INTERP)))
    for i in range(irred2_ifcs.shape[1]):
        new_irred2_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred2_ifcs[:,i])

    if p.rc3 is not None:
        irred3_ifcs = np.array(irred3_ifcs)
        new_irred3_ifcs = np.zeros((irred3_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred3_ifcs.shape[1]):
            new_irred3_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred3_ifcs[:,i])

    if p.rc4 is not None:
        irred4_ifcs = np.array(irred4_ifcs)
        new_irred4_ifcs = np.zeros((irred4_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred4_ifcs.shape[1]):
            new_irred4_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred4_ifcs[:,i])

    if p.interpolate_U0:
        interpolated_U0s = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, U0s)

    # Save interpolated IFCs per tempearture
    irred2_paths = []
    irred3_paths = [] if p.rc3 is not None else None
    irred4_paths = [] if p.rc4 is not None else None
    irred_out_path = join(paths.basepath, "INTERPOLATED_IRRED_IFCS")
    os.mkdir(irred_out_path)
    for j, T in enumerate(TEMPS_INTERP):
        temp_str = temp_to_str(T)
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
                     fmt = "%.12f", header = "# Temperature U0")
    #* TODO INTERPOLATE SECOND ORDER CUMULANT

    # Writes out command that must be used to get the same forcemap
    ef_tmp = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    with open(join(paths.basepath, "READ_IRRED_IFC_CMD.txt"), "w") as f:
        f.write("""# YOU MUST USE THIS COMMAND WHEN USING THE INTERPOLATED IFCS.
                # THE PARAMS MUST MATCH THOSE USED TO CREATE THE INTERPOLATION POINTS.
                # THE infile.ucposcar AND infile.ssposcar MUST ALSO BE IDENTICAL."""
                )
        f.write(ef_tmp._cmd())


    if p.make_ss_ifcs:
        logging.info("Converting interpolated irreducible IFCs, back to super cell IFCs")
        ss_out_dir = convert_irred_to_ss_ifc(p, paths, N_atoms, TEMPS_INTERP, TEMPS_SIM, sim_dir, irred_out_path)
        shutil.copyfile(join(irred_out_path, "interpolated_U0s.txt"), join(ss_out_dir, "interpolated_U0s.txt"))

    if p.cleanup:
        os.chdir(paths.basepath)
        os.system(f"rm -rf {sim_dir}")

    return irred_out_path, ss_out_dir
    
def convert_irred_to_ss_ifc(p : InterpolateIFCParams, paths : Paths, N_atoms :int,
                             TEMPS_INTERP, TEMPS_SIM, sim_dir, irred_out_path):
    
    ss_out_dir = join(paths.basepath, "INTERPOLATED_SUPERCELL_IFCS")
    os.mkdir(ss_out_dir)

    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    for T in TEMPS_INTERP:
        temp_str = temp_to_str(T)
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
        shutil.copyfile(join(paths.basepath, "infile.ucposcar"), join(sim_dir_T, "infile.ucposcar"))
        shutil.copyfile(join(paths.basepath, "infile.ssposcar"), join(sim_dir_T, "infile.ssposcar"))

        # forces, positions and stat file will not affect
        # IFC calculation but need to be valid. Just
        # use the files from the first temp
        shutil.copyfile(join(sim_dir, f"T{TEMPS_SIM[0]}", "infile.positions"), join(sim_dir_T, "infile.positions"))
        shutil.copyfile(join(sim_dir, f"T{TEMPS_SIM[0]}", "infile.forces"), join(sim_dir_T, "infile.forces"))
        stat_fmt = ['%d', '%d'] + ['%.8f'] * 11
        np.savetxt(join(sim_dir_T, "infile.stat"), np.zeros((p.lds.n_configs, 13)), fmt = stat_fmt) # this isnt even used, but still has to exist 

        # Write a new meta file cause we can
        write_tdep_meta(sim_dir_T, N_atoms, p.lds.n_configs, p.lds.time_step_fs, T)

        # Re-run extract force constants using the irreducible IFCs
        # This still needs infile.meta, infile.stat, infile.positions and infile.forces
        # to run although those will not be used for calculating the IFCs
        ef.mpirun(p.n_cores_max, sim_dir_T)

        # Move force constants to make folder to avoid clean up 
        shutil.copyfile(join(sim_dir_T, "outfile.forceconstant"), join(ss_out_dir, f"infile.forceconstant_{temp_str}"))
        
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_thirdorder"), join(ss_out_dir, f"infile.forceconstant_thirdorder_{temp_str}"))
            
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_fourthorder"), join(ss_out_dir, f"infile.forceconstant_fourthorder_{temp_str}"))

    return ss_out_dir
