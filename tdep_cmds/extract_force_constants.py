from typing import Optional
import numpy as np
import os

from cmds.tdep_cmd import TDEP_Command

class ExtractForceConstants(TDEP_Command):
    
    def __init__(
            self,
            r_cut2 : float,
            r_cut3 : Optional[float] = None,
            r_cut4 : Optional[float] = None,
            stride : Optional[int] = 1,
            log_file : str = "tdep.log"
        ):
        self.r_cut2 = r_cut2
        self.r_cut3 = r_cut3
        self.r_cut4 = r_cut4
        self.stride = stride
        self.log_file = log_file
    
    @property
    def input_files(self):
        return ["infile.ucposcar", "infile.ssposcar", "infile.meta", "infile.stat", "infile.positions", "infile.forces"]

    def input_parameters_valid(self):

        if self.stride <= 0 or not isinstance(self.stride, int):
            raise ValueError(f"stride in extract_forceconstants must be an integer greater than 0. Got : {self.stride}.")
        
        rcs = {"rc2" : self.r_cut2, "rc3" : self.r_cut3, "rc4" : self.r_cut4}
        if np.any([rc for rc in rcs.values() if rc is not None] < 0.0):
            raise ValueError(f"Force constant cutoffs in extract_forceconstants must be positive. Got : {rcs}")
        
        return True

    def __cmd(self) -> str:
        cmd = f"extract_forceconstants -s {self.stride} -rc2 {self.r_cut2}"

        if self.r_cut3 is not None:
            cmd += f" -rc3 {self.r_cut3}"
        
        if self.r_cut4 is not None:
            cmd += f" -rc4 {self.r_cut4}"

        cmd += f"--potential_energy_differences --verbose > {self.log_file}"

        return cmd