from dataclasses import dataclass
from typing import Optional
import logging
import os, shutil
from os.path import join

from src import (
    PathLike,
    LammpsCalculator,
    CanonicalConfigs,
    ExtractForceConstants, 
    write_tdep_meta
)

@dataclass
class sTDEP_Params:
    iters : int
    n_configs : int
    temperature : float
    mode : str
    basepath : PathLike 
    maximum_frequency : float
    r_cut2 : float
    ncores : int = 1
    force_calc : str = "lammps"
    lammps_base_script : Optional[str] = None #(e.g. LJ_argon.in)

def run_init_iteration(p : sTDEP_Params, current_dir : PathLike, run_dir : PathLike):

    os.makedirs(run_dir)
    shutil.copyfile(join(current_dir, "infile.ssposcar"), join(run_dir, "infile.ssposcar"))
    shutil.copyfile(join(current_dir, "infile.ucposcar"), join(run_dir, "infile.ucposcar"))

    cc = CanonicalConfigs(p.mode, p.n_configs, p.temperature, maximum_frequency = p.maximum_frequency)
    res = cc.run(run_dir)

    if res != 0:
        logging.error("Possible error in first iteration of sTDEP")

def prepare_next_dir(current_dir, dest_dir, init_pass : bool = False):
    os.mkdir(dest_dir)
    shutil.copyfile(join(current_dir, "infile.ssposcar"), join(dest_dir, "infile.ssposcar"))
    shutil.copyfile(join(current_dir, "infile.ucposcar"), join(dest_dir, "infile.ucposcar"))
    if not init_pass:
        shutil.copyfile(join(current_dir, "outfile.forceconstant"), join(dest_dir, "infile.forceconstant"))
    else:
        shutil.copyfile(join(current_dir, "outfile.fakeforceconstant"), join(dest_dir, "infile.forceconstant"))

def run_stdep(p : sTDEP_Params):

    iter_path = lambda i : join(p.basepath, f"stdep_iter_{i}")
    
    # Generate force constants from maximum frequency
    # These will seed the self-consistent iteration
    init_iter_path = iter_path("init")
    run_init_iteration(p, p. basepath, iter_path("init"))

    # Seed first iteration
    prepare_next_dir(init_iter_path, iter_path(0), True)

    # Pre-calc stuff needed for force calculators
    if p.force_calc == "lammps":
        file_path = os.path.dirname(os.path.realpath(__file__))
        base_infile_path = os.path.join(file_path, "..", "..", "lammps_scripts", p.lammps_base_script)
        base_infile_path = os.path.abspath(base_infile_path)
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")

    cc = CanonicalConfigs(p.mode, p.n_configs, p.temperature)
    for i in range(p.iters):
        ip = iter_path(i)
        
        # Generate configurations
        res = cc.run(ip)

        if res != 0:
            logging.error(f"Possible error in iteration {i} of sTDEP")

        # Calculate forces on configurations
        if p.force_calc == "lammps":
            N_atoms = cc.output_to_lammps_input(ip)
            lc = LammpsCalculator(base_infile_path,
                                  ip,
                                  "contcar_conf",
                                  p.n_configs)
            lc.run(ip)
            lc.remove_dump_headers(ip) # makes infile.positions, infile.stat, infile.forces
            write_tdep_meta(ip, N_atoms, p.n_configs, 0.0, p.temperature)
        else:
            raise NotImplementedError("VASP Force calc not implemented")


        ef = ExtractForceConstants(p.r_cut2)
        ef.mpirun(p.ncores, ip)


        # make_dos_plot()
        if i < p.iters - 1:
            prepare_next_dir(ip, iter_path(i+1))

    # Using Last Dataset fit third, fourth order IFCs

