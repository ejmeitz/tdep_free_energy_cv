import numpy as np
import logging
import os, shutil
from os.path import join
from typing import List, Optional
from dataclasses import dataclass

from src import (
    ExtractForceConstants,
    LammpsSimulator,
    PathLike,
    setup_logging,
    remove_dump_headers,
    write_tdep_meta
)

@dataclass
class LammpsDynamicsSettings:
    init_structure : str #(e.g. initial_structure_LJ.data)
    time_step : float 
    n_atoms : int
    n_cores : int = 4
    data_interval : int = 5000
    n_configs : int = 500
    thermo_data_temp_idx : int = 1

@dataclass
class InterpolateIFCParams:
    basepath : PathLike
    temps_to_calculate : List[float]
    temps_to_interpolate : List[float]
    rc2 : float
    rc3 : Optional[float] = None
    rc4 : Optional[float] = None
    n_cores : int = 1
    # interp_mode : str = "linear" #only piecewise linear for now
    force_calc : str = "lammps"
    lammps_base_script : Optional[str] = None #(e.g., LJ_argon_dynamics.in)
    lammps_init_structure : Optional[str] = None #(e.g. initial_structure_LJ.data)
    lds : Optional[LammpsDynamicsSettings] = None
    

def run_lammps(p, file_path, T, base_infile_path, sim_root_dir):
    structure_path = os.path.join(file_path, "..", "..", "data", "structures", "LJ", p.lammps_init_structure)
    N_steps = p.lds.n_configs * p.lds.data_interval
    var_dict = {"T" : T, "N_steps" : N_steps, "structure_path" : structure_path}
    ls = LammpsSimulator(base_infile_path, sim_root_dir, var_dict)

    ls.mpirun(sim_root_dir, p.lds.n_cores)
    remove_dump_headers(sim_root_dir) # makes infile.positions, infile.stat, infile.forces
    write_tdep_meta(sim_root_dir, p.lds.n_atoms, p.lds.n_configs, p.lds.time_step, T)

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

def run_interpolate_ifcs(p : InterpolateIFCParams):

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
        base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.lammps_base_script)
        base_infile_path = os.path.abspath(base_infile_path)
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")

    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4)

    sim_dir = os.mkdir(join(p.basepath, "simulation_data"))

    mean_sim_temps = []
    irred2_ifcs = []
    irred3_ifcs = []
    irred4_ifcs = []
    for i, T in enumerate(TEMPS_CALC):
        sim_root_dir = join(sim_dir, f"T{T}")
        os.mkdir(sim_root_dir)

        if p.force_calc == "lammps":
            mean_sim_temp = run_lammps(p, file_path, T, base_infile_path, sim_root_dir)
            mean_sim_temps.append(mean_sim_temp)
        else:
            raise NotImplementedError()

        # Calculate force constants from lammps simulation
        res = ef.mpirun(p.n_cores, sim_root_dir)
        if res != 0:
            logging.error(f"Possible error when extracting IFCs")

        # Store irred ifcs
        tmp = np.loadtxt(join(sim_root_dir, f"infile.irrifc_secondorder"))
        irred2_ifcs.append(tmp)

        if p.rc3 is not None:
            tmp = np.loadtxt(join(sim_root_dir, f"infile.irrifc_thirdorder"))
            irred3_ifcs.append(tmp)
        
        if p.rc4 is not None:
            tmp = np.loadtxt(join(sim_root_dir, f"infile.irrifc_fourthorder"))
            irred4_ifcs.append(tmp)

    # Save actual temps to see if our IFCs are accurate
    if p.force_calc == "lammps":       
        write_temps_file(p.basepath, TEMPS_CALC, mean_sim_temps)

    # Interpolate the irreducible IFCs
    irred2_ifcs = np.array(irred2_ifcs)
    new_irred2_ifcs = np.zeros((irred2_ifcs.shape[1], len(TEMPS_INTERP)))
    for i in range(irred2_ifcs.shape[1]):
        new_irred2_ifcs[i,:] = np.interp(TEMPS_INTERP, TEMPS_CALC, irred2_ifcs[:,i])

    if p.rc3 is not None:
        irred3_ifcs = np.array(irred3_ifcs)
        new_irred3_ifcs = np.zeros((irred3_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred2_ifcs.shape[1]):
            new_irred3_ifcs[i,:] = np.interp(TEMPS_INTERP, TEMPS_CALC, irred3_ifcs[:,i])

    if p.rc4 is not None:
        irred4_ifcs = np.array(irred4_ifcs)
        new_irred4_ifcs = np.zeros((irred4_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred2_ifcs.shape[1]):
            new_irred4_ifcs[i,:] = np.interp(TEMPS_INTERP, TEMPS_CALC, irred4_ifcs[:,i])

    # Save interpolated IFCs per tempearture
    for j, T in enumerate(TEMPS_INTERP):
        np.savetxt(join(p.basepath, f"infile.irrifc_secondorder_{T.replace(".", "_")}"), new_irred2_ifcs[:,j])

        if p.rc3 is not None:
            np.savetxt(join(p.basepath, f"infile.irrifc_thirdorder_{T.replace(".", "_")}"), new_irred3_ifcs[:,j])

        if p.rc4 is not None:
            np.savetxt(join(p.basepath, f"infile.irrifc_fourthorder_{T.replace(".", "_")}"), new_irred4_ifcs[:,j])
    


    
