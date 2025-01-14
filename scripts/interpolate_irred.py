import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    InterpolateIFCParams,
    run_interpolate_irred
)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing configuration. See keys in \
                            /algorithms/interpolate_ifcs InterpolateIFCParams")
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.config) as f:
        config = InterpolateIFCParams(**yaml.safe_load(f))

    if not os.path.isdir(args.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(args.basepath)

    setup_logging("sTDEP.log", args.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    run_interpolate_irred(config)


if __name__ == "__main__":
    main()