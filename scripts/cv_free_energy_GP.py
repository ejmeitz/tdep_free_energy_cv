import argparse
import logging
import numpy as np
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
    HeatCapFreeEnergyParams,
    run_cv_free_energy_GP
)

desc = \
"""
Calculates heat capacity from the second derivatve of free energy. Force constants are interpolated
to each temperature on the finite difference stencil. By default the interpolation is done at each call
but prior data can be used by setting `Paths.interp_data_path` to a directory with the output
from interpolate_irred.

Set ifc_path in the Paths configuration to use existing force constants as interpolation nodes. Expects
format produced by ifc_from_MD command. You MUST use the same force constant cutoffs, same unit cell and same supercell
as used when calcuating the origianl set of IFCs.
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

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")
    
    if paths.ifc_path is not None:
        assert os.path.isdir(paths.ifc_path)
    else:
        raise RuntimeError("Expected ifc_path pointing to IFCs created by ifc_from_MD")

    os.chdir(paths.basepath)

    setup_logging("cv_free_energy.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_cv_free_energy_GP(params, paths)

if __name__ == "__main__":
    main()