import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
    HeatCapFreeEnergyParams,
    run_cv_free_energy,
    InterpolateIFCParams,
    LammpsDynamicsSettings
)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing configuration. See keys in \
                            /algorithms/cv_free_energy HeatCapFreeEnergyParams. Must also contain ucposrcar_path and ssposcar_path.")
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

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("sTDEP.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_cv_free_energy(params, paths)

if __name__ == "__main__":
    main()