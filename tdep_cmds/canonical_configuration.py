from typing import Optional
import os
import numpy as np

from cmds.tdep_cmd import TDEP_Command

class GenerateCanonicalConfigs(TDEP_Command):

    def __init__(
            self,
            mode : str,
            n_samples : int,
            temperature : float,
            maximum_frequency : Optional[float] = None,
            log_file : str = "tdep.log"
        ):
        self.mode = mode.lower()
        self.n_samples = n_samples
        self.tempearture = temperature
        self.maximum_frequency = maximum_frequency
        self.log_file = log_file

    @property
    def input_files(self):
        if self.maximum_frequency is None:
            return ["infile.ucposcar", "infile.ssposcar", "infile.forceconstant"]
        else:
            return ["infile.ucposcar", "infile.ssposcar"]
    
    def input_parameters_valid(self):
        if self.mode not in ["quantum", "classical"]:
            raise ValueError(f"Mode in canonical_configuration must be quantum or classical. Got : {self.mode}.")
        
        if self.n_samples <= 0 or not isinstance(self.n_samples, int):
            raise ValueError(f"n_samples in canonical_configuration must be an integer greater than 0. Got : {self.n_samples}.")
        
        if self.temperature <= 0:
           raise ValueError(f"temperature in canonical_configuration must be greater than 0. Got : {self.temperature}.") 
        
        return True

    def __cmd(self) -> str:
        if self.maximum_frequency is None:
            return f"canonical_configuration -n {self.n_samples} -t {self.tempearture} > {self.log_file}"
        else:
            return f"canonical_configuration -n {self.n_samples} -t {self.tempearture} -mf {self.maximum_frequency} > {self.log_file}"
    
    def parse_output(self):
        # `fmt` was left default so this is always VASP POSCAR
        # Return np arrays of positions and velocities
        pass

    def output_to_lammps_input():
        pass

