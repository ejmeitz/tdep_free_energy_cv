from dataclasses import dataclass, field
from typing import List, Optional
import os
from ..util import PathLike


@dataclass
class Paths:
    basepath : PathLike
    ucposcar_path : PathLike
    ssposcar_path : PathLike
    ifc_path : Optional[PathLike] = None

@dataclass
class LammpsDynamicsSettings:
    time_step_fs : float 
    infile_path : str #(e.g., /home/emeitz/LJ_argon_dynamics.in)
    structure_path : str #(e.g. /home/emeitz/initial_structure_LJ.data)
    n_cores : int = 4 # cores per MD simulation
    data_interval : int = 5000
    n_configs : int = 500
    thermo_data_temp_idx : int = 1

@dataclass
class IFC_MD_Params:
    temperatures : List[float]
    n_cores_max : int # total number of cores available
    lds : LammpsDynamicsSettings
    rc2 : float
    rc3 : Optional[float] = None
    rc4 : Optional[float] = None
    cleanup : bool = True # deletes lammps simulations at end
    

@dataclass
class sTDEP_Params:
    iters : int
    n_configs : int
    temperature : float
    mode : str
    maximum_frequency : float
    rc2 : float
    n_cores : int = 1
    force_calc : str = "lammps"
    lammps_base_script : Optional[str] = None #(e.g. LJ_argon_snapshots.in)

@dataclass
class InterpolateIFCParams:
    temps_to_interpolate : List[float]
    n_cores_max : int # total number of cores available
    rc2 : float
    rc3 : Optional[float] = None
    rc4 : Optional[float] = None
    simulated_temps : Optional[List[float]] = None # only optional if ifc_path is set in Paths
    n_cores_per_sim : int = 1 # number of cores used per simulation when calculating forces
    interp_mode = "gpr" #gpr or cubic or linear
    interpolate_U0 : bool = True
    cleanup : bool = True
    force_calc : str = "lammps"
    make_ss_ifcs : bool = False # if irred are converted back to ifc for supercell
    lds : Optional[LammpsDynamicsSettings] = None
    

@dataclass
class HeatCapFreeEnergyParams:
    temperatures : List[float]
    dT : float
    k_mesh : List[int]
    interp_settings : InterpolateIFCParams = field(
        default_factory=InterpolateIFCParams
    )
    n_cores : int = 1
    fd_stencil_size : int = 3
    quantum : bool = False # use Bose-Einstein or Classical Occupation
    stochastic : bool = False # if 2nd order IFCs are fit with sTDEP then True
    force_calc : str = "lammps"
