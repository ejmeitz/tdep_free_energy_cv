from dataclasses import dataclass
from typing import Optional, List
import logging
import os, shutil
import numpy as np
from os.path import join

from src import (
    PathLike,
    ExtractForceConstants, 
    InterpolateIFCParams,
    write_tdep_meta,
    run_interpolate_irred
)

@dataclass
class HeatCapFreeEnergyParams:
    temperature : float
    dT : float
    k_mesh : List[int]
    interp_settings : InterpolateIFCParams
    fd_stencil_size : int = 3

# Heat capacity is second derivative of free energy
def get_fd_coeffs(size):
    if size == 3:
        return [1, -2, 1], [-1, 0, 1]
    elif size == 5:
        return [-1/12, 4/3, -5/2, 4/3, -1/12], [-2, -1, 0, 1, 2]
    elif size == 7:
        return [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90], [-3,-2,-1,0,1,2,3]
    else:
        raise ValueError(f"Only support stencils of size 3,5,7 for second-order central diff. Got {size}")

def check_params_consistent(p : HeatCapFreeEnergyParams, temp_stencil : List[float]) -> None:

    if not set(temp_stencil).issubset(set(p.interp_settings.temps_to_interpolate)):
        raise ValueError("Interpolation settings inconsistent with FD stencil. Not all temps on stencil calculated by interpolation.")

    if p.interp_settings.convert_ifcs == False:
        logging.warning("Setting `convert_ifcs` to True in interp_settings. This is required for anharmonic_free_energy calculation.")
        p.interp_settings.convert_ifcs = True

def run_cv_free_energy(p : HeatCapFreeEnergyParams):


    coeffs, temp_offsets = get_fd_coeffs(p.fd_stencil_size)
    temp_stencil = p.temperature + np.array([to*p.dT for to in temp_offsets])

    check_params_consistent(p, temp_stencil)

    # Run interpolation, expects infile.ucposcar and infile.ssposcar in root dir
    run_interpolate_irred(p.interp_settings)

