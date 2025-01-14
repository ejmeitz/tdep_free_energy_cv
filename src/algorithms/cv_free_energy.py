from dataclasses import dataclass
from typing import Optional, List
import logging
import os, shutil
from os.path import join

from src import (
    PathLike,
    LammpsSimulator,
    CanonicalConfigs,
    ExtractForceConstants, 
    PhononDispersion,
    write_tdep_meta,
    remove_dump_headers
)

@dataclass
class HeatCapFreeEnergyParams:
    temperature : float
    fd_order : int
    dT : float
    k_mesh : List[int]
