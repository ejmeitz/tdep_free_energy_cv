import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    InterpolateIFCParams,
    LammpsDynamicsSettings,
    run_interpolate_irred
)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing configuration. See keys in \
                            /algorithms/interpolate_ifcs InterpolateIFCParams. Must also contain ucposrcar_path and ssposcar_path.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)
        interp_params  = {k:v for k,v in cfg_data.items() if k in InterpolateIFCParams.__annotations__}
        config = InterpolateIFCParams(**interp_params)
        if config.lds is not None:
            config.lds = LammpsDynamicsSettings(**config.lds)

    if not os.path.isdir(config.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(config.basepath)

    setup_logging("sTDEP.log", config.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(cfg_data["ucposcar_path"], os.path.join(config.basepath, "infile.ucposcar"))
    shutil.copyfile(cfg_data["ssposcar_path"], os.path.join(config.basepath, "infile.ssposcar"))

    run_interpolate_irred(config)


if __name__ == "__main__":
    main()