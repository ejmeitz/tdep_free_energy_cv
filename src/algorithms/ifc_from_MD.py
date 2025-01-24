import numpy as np
import logging
import os, shutil
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

def run_lammps(p, T, infile_path, sim_root_dir):
    N_steps = p.lds.n_configs * p.lds.data_interval
    var_dict = {"T" : T, "N_steps" : N_steps, "structure_path" : p.lds.structure_path}
    required_vars = ["T", "N_steps", "structure_path"]
    ls = LammpsSimulator(infile_path, sim_root_dir, var_dict, required_vars)

    ls.mpirun(sim_root_dir, p.lds.n_cores)
    N_atoms = get_n_atoms_from_dump(join(sim_root_dir, "dump.positions"))

    remove_dump_headers(sim_root_dir) # makes infile.positions, infile.stat, infile.forces
    write_tdep_meta(sim_root_dir, N_atoms, p.lds.n_configs, p.lds.time_step_fs, T)

    # Parse thermodata.txt and get actual average temp
    temps = np.loadtxt(join(sim_root_dir, "thermo_data.txt"))[: , p.lds.thermo_data_temp_idx]
    return np.mean(temps), N_atoms

def write_temps_file(out_dir, mean_sim_temps : dict):

    temps = np.sort(list(mean_sim_temps.keys()))

    actual_temps_path = join(out_dir, "ACTUAL_SIM_TEMPS.txt")
    with open(actual_temps_path, "w") as f:
        f.write("# Requested Temperatures:\n")
        for T in temps:
            f.write(f"{T} ")
        f.write("\n")
        f.write("# Mean Simulation Temperatures:\n")
        for T in temps:
            f.write(f"{mean_sim_temps[T]} ")
        f.write("\n")

    return actual_temps_path

def work_unit(T, *, p, infile_path, sims_dir):
    temp_str = temp_to_str(T)
    sim_root_dir = join(sims_dir, f"T{temp_str}")
    os.mkdir(sim_root_dir)
    logging.info(f"Running LAMMPS at {T} K")
    mean_sim_temp, N_atoms = run_lammps(p, T, infile_path, sim_root_dir)
    return mean_sim_temp, T

def run_ifc_from_MD(p : IFC_MD_Params, paths : Paths) -> None:
    
    sim_dir = join(paths.basepath, "simulation_data")
    os.mkdir(sim_dir)

    ifc_dir = join(paths.basepath, "IFCs")
    os.mkdir(ifc_dir)

    mean_sim_temps = {}
    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4)

    # Run LAMMPS first in parallel
    max_lammps_procs = int(np.floor(p.n_cores_max / p.lds.n_cores))
    work = partial(work_unit, p = p, infile_path = p.lds.infile_path, sims_dir = sim_dir)
    with ProcessPoolExecutor(max_workers = max_lammps_procs) as exec:
        futures = [exec.submit(work, T) for T in p.temperatures]

        for f in as_completed(futures):
            sim_temp, actual_temp = f.result()
            logging.info(f"LAMMPS Complete for Temperature: {actual_temp}K")
            mean_sim_temps[sim_temp] = actual_temp

    # Go back into each dir and extract IFCs
    for T in p.temperatures:
        temp_str = temp_to_str(T)
        sim_root_dir_T = join(sim_dir, f"T{temp_str}")

        shutil.copyfile(join(paths.basepath, "infile.ucposcar"), join(sim_root_dir_T, "infile.ucposcar"))
        shutil.copyfile(join(paths.basepath, "infile.ssposcar"), join(sim_root_dir_T, "infile.ssposcar"))

        res = ef.mpirun(p.n_cores_max, sim_root_dir_T)
        if res != 0:
            logging.error(f"Possible error when extracting IFCs")

        ifc_dir_T = join(ifc_dir, f"T{temp_str}")
        os.mkdir(ifc_dir_T)
        # COPY IFCS to output dir
        shutil.copyfile(join(sim_root_dir_T, "outfile.forceconstant"),
                        join(ifc_dir_T, f"infile.forceconstasnt"))
        shutil.copyfile(join(sim_root_dir_T, "outfile.irrifc_secondorder"),
                        join(ifc_dir_T, f"infile.irrifc_secondorder"))
        if p.rc3 is not None:
            shutil.copyfile(join(sim_root_dir_T, "outfile.forceconstant_thirdorder"),
                        join(ifc_dir_T, f"infile.forceconstasnt_thirdorder"))
            shutil.copyfile(join(sim_root_dir_T, "outfile.irrifc_thirdorder"),
                        join(ifc_dir_T, f"infile.irrifc_thirdorder"))
        if p.rc4 is not None:
            shutil.copyfile(join(sim_root_dir_T, "outfile.forceconstant_fourthorder"),
                        join(ifc_dir_T, f"infile.forceconstasnt_fourthorder"))
            shutil.copyfile(join(sim_root_dir_T, "outfile.irrifc_fourthorder"),
                        join(ifc_dir_T, f"infile.irrifc_fourthorder"))
            
        shutil.copyfile(join(sim_root_dir_T, ef.log_file), join(ifc_dir_T, ef.log_file))
        shutil.copyfile(join(sim_root_dir_T, "outfile.U0"), join(ifc_dir_T, "outfile.U0"))
        shutil.copyfile(join(sim_root_dir_T, "outfile.energies"), join(ifc_dir_T, "outfile.energies"))

    # Only deletes useless duplicate files
    # TDEP infiles have the same info 
    # and are useful later    
    if p.cleanup:
        os.chdir(sim_dir)
        os.system(f"rm -rf dump.*")

    write_temps_file(paths.basepath,  mean_sim_temps)


    return sim_dir, ifc_dir, mean_sim_temps