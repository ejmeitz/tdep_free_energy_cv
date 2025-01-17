import logging
import os, shutil
from os.path import join

from src import (
    PathLike,
    LammpsSimulator,
    CanonicalConfigs,
    ExtractForceConstants, 
    PhononDispersion,
    write_tdep_meta,
    remove_dump_headers
)

from .configs import sTDEP_Params

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

    iter_path = lambda i : join(p.basepath, f"stdep_iter_{i}") if i >= 0 else join(p.basepath, "stdep_iter_init")
    
    # Generate force constants from maximum frequency
    # These will seed the self-consistent iteration
    init_iter_path = iter_path(-1)
    run_init_iteration(p, p.basepath, init_iter_path)

    # Seed first iteration
    prepare_next_dir(init_iter_path, iter_path(0), True)

    # Pre-calc stuff needed for force calculators
    if p.force_calc == "lammps":
        file_path = os.path.dirname(os.path.realpath(__file__))
        base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.lammps_base_script)
        base_infile_path = os.path.abspath(base_infile_path)
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")

    ## INITIALIZE TDEP COMMANDS ##

    cc = CanonicalConfigs(p.mode, p.n_configs, p.temperature)

    # 40x40x40 k-grid should be more than enough for any material
    # I dont really feel like parsing a list as a cmd argument
    pd = PhononDispersion([40, 40, 40])

    ef = ExtractForceConstants(p.r_cut2)
    
    for i in range(p.iters):
        ip = iter_path(i)
        
        # Generate configurations
        res = cc.run(ip)

        if res != 0:
            logging.error(f"Possible error in iteration {i} of sTDEP")

        # Calculate forces on configurations
        if p.force_calc == "lammps":
            N_atoms = cc.output_to_lammps_input(ip)

            # These must match the keys in the LAMMPS input file
            var_dict = {"structure_path" : os.path.join(ip, "contcar_conf"), "n_configs" : p.n_configs}
            ls = LammpsSimulator(base_infile_path, ip, var_dict)

            ls.run(ip)
            remove_dump_headers(ip) # makes infile.positions, infile.stat, infile.forces
            write_tdep_meta(ip, N_atoms, p.n_configs, 0.0, p.temperature)
        else:
            raise NotImplementedError("VASP Force calc not implemented")

        # Calculate force constants
        res = ef.mpirun(p.ncores, ip)

        if res != 0:
            logging.error(f"Possible error in iteration {i} of sTDEP when extracting IFCs")

         # Calculate Dispersion and DOS
        res = pd.run(ip)

        if res != 0:
            logging.error(f"Possible error in iteration {i} of sTDEP when calculating dispersion")

        # infile.forceconstant is actually from the last iteration
        # set outpath to be previous iteration
        pd.plot_dos(ip, iter_path(i-1))
        pd.plot_dispersion(ip, iter_path(i-1))

        if i < p.iters - 1:
            prepare_next_dir(ip, iter_path(i+1))

    # Make dir for final results
    results_path = os.path.join(p.basepath, "RESULTS")
    os.mkdir(results_path)
    # Copy final IFCs to root dir
    final_iter_path = iter_path(p.iters - 1)
    shutil.copyfile(join(final_iter_path, "outfile.forceconstant"),
                    join(results_path, "infile.forceconstant"))
    shutil.copyfile(join(p.basepath, "infile.ucposcar"),
                    join(results_path, "infile.ucposcar"))

    # Plot final DOS and Dispersion
    res = pd.run(results_path)
    if res != 0:
        logging.error(f"Possible error calculating dispersion for last iteration of sTDEP")
    pd.plot_dos(results_path, results_path)
    pd.plot_dispersion(results_path, results_path)
    # Also save in correct folder to avoid confusion
    pd.plot_dos(results_path, final_iter_path)
    pd.plot_dispersion(results_path, final_iter_path)

    # Using Last Dataset fit third, fourth order IFCs????

    # Make dir for convergence studies
    conv_path = os.path.join(p.basepath, "CONVERGENCE")
    os.mkdir(conv_path)
    # Make DOS convergence plot
    # Make R^2 convergence plot
    

