import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
    run_ifc_from_MD,
    LammpsDynamicsSettings,
    IFC_MD_Params
)

"""
Calculates force constants at a pre-determined set of temperatures. Automatically runs LAMMPS
to generate configurations.
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing configuration.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)

    paths = Paths(**cfg_data["Paths"])
    params = IFC_MD_Params(**cfg_data["IFC_MD_Params"])

    if params.interp_settings.lds is not None:
        params.interp_settings.lds = LammpsDynamicsSettings(**params.interp_settings.lds)
    else:
        raise ValueError("Expected key lds in IFC_MD_Params")

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("sTDEP.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_ifc_from_MD(params, paths)


if __name__ == "__main__":
    main()