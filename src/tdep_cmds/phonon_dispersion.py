from typing import List
import matplotlib.pyplot as plt
import os

from .tdep_cmd import TDEP_Command
from src import PathLike

class PhononDispersion(TDEP_Command):

    def __init__(
            self,
            q_grid : List[int],
            dos : bool = True,
            gruneisen : bool = False,
            unit = "thz",
            log_file : str = "phonon_disp.log"
        ):
        self.q_grid = q_grid
        self.dos = dos
        self.gruneisen = gruneisen
        self.unit = unit
        self.log_file = log_file


    @property
    def input_files(self):
        if self.gruneisen:
            return ["infile.ucposcar", "infile.forceconstant", "infile.forceconstant_thirdorder"]
        else:
            return ["infile.ucposcar", "infile.forceconstant"]
    
    def output_files(self):
        files = ["outfile.dispersion_relations.hdf5", "outfile.dispersion_relations", "oufile.group_velocities"]
        if self.dos:
            files.append("outfile.phonon_dos.hdf5")
            files.append("outfile.phonon_dos")
        if self.gruneisen:
            files.append("outfile.mode_gruneisen_parameters")
    
    def input_parameters_valid(self):
        if not all([isinstance(q, int) for q in self.q_grid]):
            raise ValueError(f"All q-grid entries must be integers. Got : {self.q_grid}.")
        
        if len(self.q_grid) != 3:
            raise ValueError(f"q-grid must have length 3. Got : {len(self.q_grid)}")
        
        if self.unit not in ["thz", "mev", "icm"]:
            raise ValueError(f"Unit of in phonon_dispersion must be thz, mev or icm. Got : {self.unit}")
        
        return True

    def _cmd(self) -> str:
        cmd = f"phonon_dispersion_relations -qg {self.q_grid[0]} {self.q_grid[1]} {self.q_grid[2]}"
        if self.dos:
            cmd += " --dos"
        
        if self.gruneisen:
            cmd += " --gruneisen"

        return cmd + f" > {self.log_file}"

    def plot_dos(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        pass

    def plot_dispersion(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        pass

    def plot_gruneisen(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        pass