import argparse
import logging
import numpy as np
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
    HeatCapFreeEnergyParams,
    run_cv_free_energy,
    initialize_free_energy, 
    run_interpolate_irred,
    InterpolateIFCParams,
    LammpsDynamicsSettings,
    simulate_interp_node_data
)

desc = \
"""
Calculates heat capacity from the second derivatve of free energy. Force constants are interpolated
to each temperature on the finite difference stencil. By default the interpolation is done at each call
but prior data can be used by setting `Paths.interp_data_path` to a directory with the output
from interpolate_irred.

Set ifc_path in the Paths configuration to use existing force constants as interpolation nodes. Expects
format produces by ifc_from_MD command.
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing Paths and InterpolateIFCParams configuration.")
    return parser.parse_args()

def main():

    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)

    paths = Paths(**cfg_data["Paths"])
    params = HeatCapFreeEnergyParams(**cfg_data["HeatCapFreeEnergyParams"])
    if params.interp_settings is not None:
        params.interp_settings = InterpolateIFCParams(**params.interp_settings)
    else:
        raise ValueError("Expected key interp_settings in HeatCapFreeEnergyParams")

    if params.interp_settings.lds is not None:
        params.interp_settings.lds = LammpsDynamicsSettings(**params.interp_settings.lds)
    else:
        if paths.ifc_path is not None:
            raise ValueError("No LAMMPS settings passed, expected ifc_path to be set in Paths section.")

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("cv_free_energy.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_interpolate = paths.ifc_path is None
    temp_stencil = initialize_free_energy(params, run_interpolate)

    # Run simulations at interpolation nodes if necessary
    sim_dir = None; ifc_dir = None; mean_sim_temps = None
    if paths.ifc_path is not None:
        sim_dir, ifc_dir, mean_sim_temps = simulate_interp_node_data(params, paths)
    else:
        # Parse existing simulation data
        sim_dir = os.path.join(paths.basepath, "simulation_data")
        os.mkdir(sim_dir)
        ifc_dir = paths.ifc_path
        mean_sim_temps = np.loadtxt(os.path.join(paths.ifc_path, "..", "ACTUAL_SIM_TEMPS.txt"))
    
    _, ss_out_dir = run_interpolate_irred(params, paths, sim_dir, ifc_dir, mean_sim_temps)


    run_cv_free_energy(params, paths, temp_stencil, ss_out_dir)

if __name__ == "__main__":
    main()