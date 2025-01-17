import argparse
import logging
import os, shutil
from os.path import join
import yaml
from src import sTDEP_Params, Paths, run_stdep, setup_logging


def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing sTDEP params. See keys in \
                            /algorithms/configs.py sTDEP_Params. Must also entry for Paths.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)
    
    params = sTDEP_Params(**cfg_data["sTDEP_Params"])
    paths = Paths(**cfg_data["Paths"])

    if not os.path.isdir(args.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(args.basepath)

    setup_logging("sTDEP.log", args.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    shutil.copyfile(paths.ucposcar_path, join(paths.basepath, "infile.ssposcar"))
    shutil.copyfile(paths.ssposcar_path, join(paths.basepath, "infile.ssposcar"))       


    run_stdep(params, paths)

if __name__ == "__main__":
    main()



# sTDEP --iters 2 --n_configs 200 --temperature 50 --mode "classical" --maximum_frequency 2.0 --basepath "/home/emeitz/tests/tdep_outputs/sTDEP_TEST" --r_cut2 8.5 --lammps_base_script "LJ_argon.in" --ncores 4