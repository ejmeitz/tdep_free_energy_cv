import argparse
import logging
import os, shutil
import yaml
import numpy as np

from src import (
    setup_logging,
    Paths,
    InterpolateIFCParams,
    LammpsDynamicsSettings,
    run_interpolate_irred,
    simulate_interp_node_data
)

desc =\
"""
Interpolates the irreducible force constants with respect to temperature. Can calculate new force constants to use
as interpolation nodes or load previously calculated force constants from the ifc_from_MD script.

To load previously calculated force constants set ifc_path in the Paths section of the config.yml file to the IFCs directory
created by ifc_from_MD.
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing Paths and InterpolateIFCParams configurations.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)

    paths = Paths(**cfg_data["Paths"])

    params = InterpolateIFCParams(**cfg_data["InterpolateIFCParams"])
    if params.lds is not None:
        params.lds = LammpsDynamicsSettings(**params.lds)
    else:
        if paths.ifc_path is not None:
            raise ValueError("No LAMMPS settings passed, expected ifc_path to be set in Paths section.")

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("IFC_from_MD.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))
    
    sim_dir = None; ifc_dir = None; mean_sim_temps = None
    if paths.ifc_path is not None:
        sim_dir, ifc_dir, mean_sim_temps = simulate_interp_node_data(params, paths)
    else:
        sim_dir = os.path.join(paths.basepath, "simulation_data")
        os.mkdir(sim_dir)
        ifc_dir = paths.ifc_path
        mean_sim_temps = np.loadtxt(os.path.join(paths.ifc_path, "..", "ACTUAL_SIM_TEMPS.txt"))
    
    run_interpolate_irred(params, paths, sim_dir, ifc_dir, mean_sim_temps)

if __name__ == "__main__":
    main()