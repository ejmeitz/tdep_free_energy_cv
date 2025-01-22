import logging
import os
from pathlib import Path
from typing import Union

PathLike = Union[Path, str]

def setup_logging(logname : str, savepath : PathLike = os.getcwd()):
    logging.basicConfig(
        level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(savepath, logname)), 
            logging.StreamHandler()
        ]
    )
    logging.info("Logger initialized")


def write_tdep_meta(out_dir : PathLike, N_atoms : int, N_samples : int, dt_fs : float, temperature : float):
    with open(os.path.join(out_dir, "infile.meta"), "w") as f:
        f.write(f"{N_atoms} # N atoms\n")
        f.write(f"{N_samples} # N samples\n")
        f.write(f"{dt_fs} # timestep in fs\n")
        f.write(f"{temperature} # temperature in K")


def temp_to_str(T: float):
    return f"{float(T)}".replace('.', '_')