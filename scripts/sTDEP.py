import argparse
import logging
import os, shutil
from src import sTDEP_Params, run_stdep, GenerateStructure, setup_logging


def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")
    
    parser.add_argument("--iters", '-i', type=int, required=True, help="Number of self-consistent iterations.")
    parser.add_argument("--n_configs", '-nc', type=int, required=True, help="Number of configurations.")
    parser.add_argument("--temperature", '-T', type=float, required=True, help="Simulation temperature.")
    parser.add_argument("--mode", type=str, required=True, default = "classical", help="Classical or Quantum.")
    parser.add_argument("--basepath", type=str, required=True, help="Base path for data storage.")
    parser.add_argument("--maximum_frequency", '-mf', type=float, required=True, help="Maximum frequency.")
    parser.add_argument("--r_cut2", type=float, required=True, help="Cutoff radius.")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores. Default is 1.")
    parser.add_argument("--force_calc", type=str, default="lammps", help="Force calculation method. Default is 'lammps'.")
    parser.add_argument("--lammps_base_script", required = False, type = str, help="Lammps base script (e.g., LJ_argon.in).")
    parser.add_argument("--num_unit_cells" , required = False, type = int, help = "If not passed expects infile.sspocar, else will generate this file")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    if not os.path.isdir(args.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(args.basepath)

    setup_logging("sTDEP.log", args.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Check if ucposcar exists
    if not os.path.isfile(os.path.join(args.basepath, "infile.ucposcar")):
        raise RuntimeError(f"Could not find ucposcar in {args.basepath}")
    
    if args.num_unit_cells is None and not os.path.isfile(os.path.join(args.basepath, "infile.ssposcar")):
        raise RuntimeError(f"num_unit_cells not passed, expected infile.ssposcar in {args.basepath}")
    elif args.num_unit_cells is not None:
        logging.info("Generating {args.num_unit_cells}x{args.num_unit_cells}x{args.num_unit_cells} supercell ")
        gc = GenerateStructure([args.num_unit_cells, args.num_unit_cells, args.num_unit_cells])
        res = gc.run(args.basepath)
        if res != 0:
            logging.error("Possible error running generate_structure")
        shutil.copyfile("outfile.ssposcar", "infile.ssposcar")        

    stdep_args = {k:v for k,v in vars(args).items() if k in sTDEP_Params.__annotations__}
    params = sTDEP_Params(**stdep_args)
    run_stdep(params)

if __name__ == "__main__":
    main()



# sTDEP --iters 2 --n_configs 200 --temperature 50 --mode "classical" --maximum_frequency 2.0 --basepath "/home/emeitz/tests/tdep_outputs/sTDEP_TEST" --r_cut2 8.5 --lammps_base_script "LJ_argon.in" --ncores 4