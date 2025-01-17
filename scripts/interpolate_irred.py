import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
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

    paths = Paths(**cfg_data["Paths"])

    params = InterpolateIFCParams(**cfg_data["InterpolateIFCParams"])
    if params.lds is not None:
        params.lds = LammpsDynamicsSettings(**params.lds)

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("sTDEP.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_interpolate_irred(params, paths)


if __name__ == "__main__":
    main()