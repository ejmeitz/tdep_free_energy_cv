from typing import Optional
import os
import glob
import logging
import numpy as np

from cmds.tdep_cmd import TDEP_Command
from src import PathLike

class CanonicalConfigs(TDEP_Command):

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
    
    def output_files(self, path):
        return glob.glob(os.path.join(path, "contcar_conf*"))
    
    def input_parameters_valid(self):
        if self.mode not in ["quantum", "classical"]:
            raise ValueError(f"Mode in canonical_configuration must be quantum or classical. Got : {self.mode}.")
        
        if self.n_samples <= 0 or not isinstance(self.n_samples, int):
            raise ValueError(f"n_samples in canonical_configuration must be an integer greater than 0. Got : {self.n_samples}.")
        
        if self.temperature <= 0:
           raise ValueError(f"temperature in canonical_configuration must be greater than 0. Got : {self.temperature}.") 
        
        if self.n_samples > 9999:
            raise ValueError(f"canonical_configuratoins is limited to 9999 configurations, tried to generate ${self.n_samples}")
        
        return True

    def __cmd(self) -> str:
        if self.maximum_frequency is None:
            return f"canonical_configuration -n {self.n_samples} -t {self.tempearture} > {self.log_file}"
        else:
            return f"canonical_configuration -n {self.n_samples} -t {self.tempearture} -mf {self.maximum_frequency} > {self.log_file}"
    
    def parse_output(self):
        # `fmt` was left default so this is always VASP POSCAR
        # Return np arrays of positions and velocities??
        # use atomsk? call fn in util.py
        pass

    def __parse_poscar_velos(self, poscar_path, out_units : str = "METAL"):

        if out_units != "METAL":
            raise ValueError("This code only supports LAMMPS metal units")

        N_atoms = None
        with open(poscar_path, "r") as f:
            for i, line in enumerate(f.readlines()):
                if i == 6:
                    N_atoms = int(line.strip())

        velocities = []
        with open(poscar_path, "r") as f:
            for i, line in enumerate(f.readlines()):
                if i == N_atoms + 9:
                        velocities.append([float(e) for e in line.strip().split()])
        
        if len(velocities) != N_atoms:
            raise RuntimeError(f"Error parsing velocities from POSCAR. Got {len(velocities)} entries but expected {N_atoms}.")

        return 10 * np.array(velocities) # POSCAR is Ang / fs, LAMMPS metal is Ang / ps

    def output_to_lammps_input(self, dir : PathLike = os.getcwd()):
        poscar_files = self.output_files(dir)
        lammps_files = []

        res = 0
        for file in poscar_files:
            res |= os.system(f"atomsk {file} lammps")
            output_file = file + ".lmp" # output of canonical config has no ext to remove
            lammps_files.append(output_file)

            # Append velocities to end of lammps file
            velos = self.__parse_poscar_velos(file)
            with open(output_file, "a") as f:
                f.write("\nVelocities\n\n")
                np.savetxt(f, velos)

        if res != 0:
            logging.error("Possible error using atomsk to convert VASP POSCAR to LAMMPS format")

        return lammps_files
