from typing import Optional
import os
import numpy as np

from tdep_cmd import TDEP_Command

class ExtractForceConstants(TDEP_Command):
    pass


class GenerateMDSamples(TDEP_Command):

    # Runs samples_from_md
    def resample_md(self):
        pass


class GenerateCanonicalConfigs(TDEP_Command):

    input_files = ["infile.ucposcar", "infile.ssposcar", "infile.forceconstant"]

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

    @property
    def input_files(self):
        return self.input_files
    
    def input_parameters_valid(self):
        if self.mode not in ["quantum", "classical"]:
            raise ValueError(f"Mode in canonical_configuration must be quantum or classical. Got : {self.mode}.")
        
        if self.n_samples <= 0:
            raise ValueError(f"n_samples in canonical_configuration must be greater than 0. Got : {self.n_samples}.")
        
        if self.temperature <= 0:
           raise ValueError(f"temperature in canonical_configuration must be greater than 0. Got : {self.temperature}.") 
        
        return True

    def __cmd(self):
        return f"canonical_configuration -n {self.n_samples} -t {self.tempearture} > {self.log_file}"
    
    def run(self, ncores : int) -> int:
        return os.system(f"mpirun -np {ncores} " + self.cmd())
    
    def parse_output(self):
        # `fmt` was left default so this is always VASP POSCAR
        # Return np arrays of positions and velocities
        pass

    def output_to_lammps_input():
        pass

