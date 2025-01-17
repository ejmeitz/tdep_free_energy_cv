import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    HeatCapFreeEnergyParams,
    run_cv_free_energy,
    InterpolateIFCParams
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
        cv_params  = {k:v for k,v in cfg_data.items() if k in HeatCapFreeEnergyParams.__annotations__}
        config = HeatCapFreeEnergyParams(**cv_params)
        if config.interp_settings is not None:
            config.interp_settings = InterpolateIFCParams(**config.interp_settings)

        if "ucposcar_path" not in cfg_data or "ssposcar_path" not in cfg_data:
            raise ValueError("Need ucposcar_path and ssposcar_path in config.yml as well")

    if not os.path.isdir(config.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(config.basepath)

    setup_logging("sTDEP.log", config.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(config.interp_settings.ucposcar_path, os.path.join(config.basepath, "infile.ucposcar"))
    shutil.copyfile(config.interp_settings.ssposcar_path, os.path.join(config.basepath, "infile.ssposcar"))

    run_cv_free_energy(config)

if __name__ == "__main__":
    main()