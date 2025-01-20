import numpy as np
import logging
import os
from os.path import join
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

from src import (
    LammpsSimulator,
    ExtractForceConstants,
    get_n_atoms_from_dump,
    remove_dump_headers,
    write_tdep_meta,
    temp_to_str
)

from .configs import IFC_MD_Params, Paths

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
    return np.mean(temps), N_atoms

def write_temps_file(out_dir, mean_sim_temps : dict):

    temps = np.sort(list(mean_sim_temps.keys()))

    with open(join(out_dir, "ACTUAL_SIM_TEMPS.txt"), "w") as f:
        f.write("# Requested Temperatures:\n")
        for T in temps:
            f.write(f"{T} ")
        f.write("\n")
        f.write("# Mean Simulation Temperatures:\n")
        for T in temps:
            f.write(f"{mean_sim_temps[T]} ")
        f.write("\n")

def work_unit(T, *, p, base_infile_path):
    temp_str = temp_to_str(T)
    sim_root_dir = join(p.basepath, f"T{temp_str}")
    logging.info(f"Running LAMMPS at {T} K")
    mean_sim_temp, N_atoms = run_lammps(p, T, base_infile_path, sim_root_dir)
    return mean_sim_temp, T

def run_ifc_from_MD(p : IFC_MD_Params, paths : Paths) -> None:
    
    file_path = os.path.dirname(os.path.realpath(__file__))
    base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.lds.script_name)
    base_infile_path = os.path.abspath(base_infile_path)

    mean_sim_temps = {}
    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4)

    # Run LAMMPS first in parallel
    max_lammps_procs = int(np.foor(p.n_cores_max / p.n_cores_per_lammps))
    work = partial(work_unit, p = p, base_infile_path = base_infile_path)
    with ProcessPoolExecutor(max_workers = max_lammps_procs) as exec:
        futures = [exec.submit(work, T) for T in p.temperatures]

        for f in as_completed(futures):
            sim_temp, actual_temp = f.result()
            logging.info(f"LAMMPS Complete for Temperature: {actual_temp}K")
            mean_sim_temps[sim_temp] = actual_temp

    # Go back into each dir and extract IFCs
    for T in p.temperatures:
        temp_str = temp_to_str(T)
        sim_root_dir = join(p.basepath, f"T{temp_str}")
        res = ef.mpirun(p.n_cores_max, sim_root_dir)
        if res != 0:
            logging.error(f"Possible error when extracting IFCs")

        if p.cleanup:
            os.chdir(sim_root_dir)
            os.system("rm -rf dump.*")

    write_temps_file(paths.basepath,  mean_sim_temps)