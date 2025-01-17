from dataclasses import dataclass, field
from typing import List, Optional
from ..util import PathLike

@dataclass
class sTDEP_Params:
    iters : int
    n_configs : int
    temperature : float
    mode : str
    basepath : PathLike 
    maximum_frequency : float
    r_cut2 : float
    ncores : int = 1
    force_calc : str = "lammps"
    lammps_base_script : Optional[str] = None #(e.g. LJ_argon_snapshots.in)


@dataclass
class LammpsDynamicsSettings:
    time_step_fs : float 
    script_name : str #(e.g., LJ_argon_dynamics.in)
    structure_path : str #(e.g. /home/emeitz/initial_structure_LJ.data)
    n_cores : int = 4
    data_interval : int = 5000
    n_configs : int = 500
    thermo_data_temp_idx : int = 1

@dataclass
class InterpolateIFCParams:
    basepath : PathLike
    temps_to_calculate : List[float]
    temps_to_interpolate : List[float]
    rc2 : float
    rc3 : Optional[float] = None
    rc4 : Optional[float] = None
    n_cores : int = 1
    interpolate_U0 : bool = True
    cleanup : bool = True
    force_calc : str = "lammps"
    make_ss_ifcs : bool = False # if irred are converted back to ifc for supercell
    lds : LammpsDynamicsSettings = field(
        default_factory=LammpsDynamicsSettings
    )
    

@dataclass
class HeatCapFreeEnergyParams:
    temperature : float
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
