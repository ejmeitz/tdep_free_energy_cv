import argparse
from src import sTDEP_Params, run_stdep


def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")
    
    parser.add_argument("--iters", type=int, required=True, help="Number of self-consistent iterations.")
    parser.add_argument("--n_configs", type=int, required=True, help="Number of configurations.")
    parser.add_argument("--temperature", type=float, required=True, help="Simulation temperature.")
    parser.add_argument("--mode", type=str, required=True, help="Classical or Quantum.")
    parser.add_argument("--basepath", type=str, required=True, help="Base path for data storage.")
    parser.add_argument("--maximum_frequency", type=float, required=True, help="Maximum frequency.")
    parser.add_argument("--r_cut2", type=float, required=True, help="Cutoff radius.")
    parser.add_argument("--ncores", type=int, default=1, help="Number of cores. Default is 1.")
    parser.add_argument("--force_calc", type=str, default="lammps", help="Force calculation method. Default is 'lammps'.")
    parser.add_argument("--lammps_base_script", required = False, type = str, help="Lammps base script (e.g., LJ_argon.in).")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    params = sTDEP_Params(**vars(args))
    run_stdep(params)

if __name__ == "__main__":
    main()