from typing import List
import logging
import numpy as np

from .tdep_cmd import TDEP_Command
from src import PathLike

class AnharmonicFreeEnergy(TDEP_Command):

    def __init__(
            self,
            q_grid : List[int],
            quantum : bool = False,
            third_order : bool = True,
            fourth_order : bool = True,
            stochastic : bool = False,
            dump_mode_values : bool = False, # requires my modified TDEP
            log_file : str = "anharm_free_energy.log"
        ):
        self.q_grid = q_grid
        self.quantum = quantum
        self.third_order = third_order
        self.fourth_order = fourth_order
        self.stochastic = stochastic
        self.dump_mode_values = dump_mode_values
        self.log_file = log_file

        if self.fourth_order and not self.third_order:
            logging.warning("In anharmonic free energy, activated third order, since fourth order was active.")
            self.third_order = True

    @property
    def input_files(self) -> List[str]:
        files = ["infile.ucposcar", "infile.ssposcar", "infile.forceconstant", \
                          "infile.stat", "infile.meta", "infile.forces", "infile.positions"]
        
        if self.third_order:
            files.append("infile.forceconstant_thirdorder")

        if self.fourth_order:
            files.append("infile.forceconstant_fourthorder")
        
        return files

    
    def output_files(self):
        return ["outfile.anharmonic_free_energy"]
    
    def input_parameters_valid(self):
        if len(self.q_grid) != 3:
            raise RuntimeError("q_grid must be length 3")
        return True
    
    # Parses the output in .log file NOT the actual output
    def parse_bulk_F_from_log(self, path : PathLike):
        data_out = []
        with open(path, "r") as f:
            all_data = f.readlines() #usually only like 1-2 kB so fine to load it all
        
        


    def _cmd(self) -> str:
        cmd =  f"anharmonic_free_energy -qg {self.q_grid[0]} {self.q_grid[1]} {self.q_grid[2]}"

        if self.third_order:
            cmd += " --thirdorder"
        if self.fourth_order:
            cmd += " --fourthorder"
        if self.quantum:
            cmd += " --quantum"
        if self.stochastic:
            cmd += " --stochastic"

        if self.dump_mode_values:
            cmd += " --modevalues"

        return cmd + f" > {self.log_file}"
    