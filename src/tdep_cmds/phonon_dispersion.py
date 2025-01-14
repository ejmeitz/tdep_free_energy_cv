from typing import List
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

from .tdep_cmd import TDEP_Command
from src import PathLike

class PhononDispersion(TDEP_Command):

    letter_map = {"GM" : "$\Gamma$", "X" : "X", "K" : "K", "R" : "R", "L" : "L"}

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
    
    def __unit_to_plot_label(self):
        if self.unit == "thz":
            return "THz"
        elif self.unit == "mev":
            return "meV"
        elif self.unit == "icm":
            return "cm^-1"
        return "unknown-unit"
    
    def __parse_path_letters_from_log(self):
        # parse for "number of paths: 4"
        # paths below that.
        pass

    def plot_dos(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        dos_data = np.loadtxt(os.path.join(basepath, "outfile.phonon_dos"))
        plt.plot(dos_data[:,0], dos_data[:,1], lw = 2.5, color = "#21deb2");
        plt.xlabel(f"Frequency [{self.__unit_to_plot_label()}]");
        plt.ylabel(f"DOS [1/{self.__unit_to_plot_label()}]");
        plt.savefig(os.path.join(outpath, "DOS.png"));
        plt.close();

    def plot_dispersion(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        filepath = os.path.join(basepath, "outfile.dispersion_relations.hdf5")
        with h5py.File(filepath, 'r') as f:
            x = f['/q_values'][:]
            xtck = f['/q_ticks'][:]
            # xtckl = f.attrs['/q_tick_labels'].decode().split()  # Decode and split into list
            y = f['/frequencies'][:]

        # Plot data
        plt.figure(1)
        plt.clf()
        plt.plot(x, y, color = "#21deb2", lw = 2)

        # Configure plot
        plt.xticks(xtck) #,xtckl)
        plt.gca().tick_params(axis='y', which='both', direction='in')  # For y minor ticks
        plt.ylabel(f"Frequency [{self.__unit_to_plot_label()}]");
        plt.xlim([0, max(x)])
        plt.savefig(os.path.join(outpath, "DISPERSION.png"));
        plt.close();

    def plot_gruneisen(self, basepath : PathLike = os.getcwd(), outpath : PathLike = os.getcwd()):
        pass