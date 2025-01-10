from dataclasses import dataclass
from typing import Optional
import logging
import os, shutil

from src import (
    PathLike,
    LammpsCalculator,
    CanonicalConfigs,
    ExtractForceConstants
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


def calc_forces_lammps():
    pass

def calc_forces_vasp():
    pass

def run_init_iteration(p : sTDEP_Params, outpath : PathLike):

    cc = CanonicalConfigs(p.mode, p.n_configs, p.temperature, maximum_frequency = p.max_freq)
    res = cc.mpirun(p.ncores, outpath)

    if res != 0:
        logging.error("Possible error in first iteration of sTDEP")

def make_dos_plot():
    pass


def prepare_next_dir(dest_dir, init_pass : bool = False):
    os.mkdir(dest_dir)
    shutil.copyfile("infile.ssposcar", os.path.join(dest_dir, "infile.ssposcar"))
    shutil.copyfile("infile.ucposcar", os.path.join(dest_dir, "infile.ucposcar"))
    if not init_pass:
        shutil.copyfile("outfile.forceconstant", os.path.join(dest_dir, "infile.forceconstant"))
    else:
        shutil.copyfile("outfile.fakeforceconstant", os.path.join(dest_dir, "infile.forceconstant"))

def run_stdep(p : sTDEP_Params):

    iter_path = lambda i : os.path.join(p.basepath, f"cc_iter_$(i)")
    
    # Generate force constants from maximum frequency
    # These will seed the self-consistent iteration
    run_init_iteration(p, iter_path("init"), use_max_freq = True)

    # Seed first iteration
    first_iter_path = iter_path(0)
    prepare_next_dir(first_iter_path)
    shutil.copyfile("outfile.fakeforceconstant", os.path.join(first_iter_path, "infile.forceconstant"))

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
        
        # Generate configurations
        res = cc.mpirun(p.ncores, iter_path(i))

        if res != 0:
            logging.error(f"Possible error in iteration {i} of sTDEP")

        # Calculate forces on configurations
        if p.force_calc == "lammps":
            cc.output_to_lammps_input(dir = iter_path(i))
            lc = LammpsCalculator(base_infile_path,
                                  iter_path(i),
                                  "contcar_conf",
                                  p.n_configs)
            lc.run()
            lc.remove_dump_headers(iter_path(i)) # makes infile.positions, infile.stat, infile.forces
        else:
            raise NotImplementedError("VASP Force calc not implemented")


        ef = ExtractForceConstants(p.r_cut2)
        ef.mpirun(p.ncores, iter_path(i))


        # make_dos_plot()
        prepare_next_dir(iter_path(i+1))

    # Using Last Dataset fit third, fourth order IFCs

