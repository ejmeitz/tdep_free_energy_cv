from abc import ABC, abstractmethod
import logging
from typing import List
from pathlib import Path
import os
import numpy as np

from src import PathLike

class TDEP_Command(ABC):

    @property
    @abstractmethod
    def input_files(self) -> List[str]:
        pass

    @abstractmethod
    def output_files(self, *args) -> List[str]:
        pass

    @abstractmethod
    def input_parameters_valid(self) -> bool:
        pass

    @property
    @abstractmethod
    def _cmd(self) -> str:
        """
        Builds command that could be called in shell. This method
        is private in favor of the `cmd` function which also runs
        several checks to alleviate debuggging TDEP issues. This
        cmd should not have `mpirun -np 4` or any other launcher.
        """
        pass

    def input_files_present(self, dir) -> bool:
        dir = Path(dir) if not isinstance(dir, Path) else dir

        if not dir.is_dir():
            raise ValueError(f"The provided path '{dir}' is not a valid directory.")
        
        missing_files = np.array([(dir / file).exists() for file in self.input_files])
        all_present = np.all(missing_files)
        if not all_present:
            logging.error(f"Missing the following files: {np.array(self.input_files)[~missing_files]}")

        return all_present
    
    def cmd(self, run_dir):
        """
        Checks if all reuqired input files are present and
        that the parameters passed to the cmd are valid. If so
        the command to run the program is generated and returned.
        """
        if self.input_files_present(run_dir) and self.input_parameters_valid():
            return self._cmd()
        else:
            raise RuntimeError("Input parameters invalid or missing proper input files. See log file for details.")
        
    def mpirun(self, ncores : int, run_dir : PathLike = os.getcwd()) -> int:
        initial_dir = os.getcwd()
        os.chdir(run_dir)
        c = f"mpirun -np {ncores} " + self.cmd(run_dir)
        logging.info(f"Running: {c} in {run_dir}")
        res = os.system(c)
        os.chdir(initial_dir)
        return res
    
    def run(self, run_dir : PathLike = os.getcwd()) -> int:
        initial_dir = os.getcwd()
        os.chdir(run_dir)
        logging.info(f"Running: {self.cmd(run_dir)} in {run_dir}")
        res = os.system(self.cmd())
        os.chdir(initial_dir)
        return res


