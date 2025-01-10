from abc import ABC, abstractmethod
from typing import List
from pathlib import Path
import os

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
    def __cmd(self) -> str:
        """
        Builds command that could be called in shell. This method
        is private in favor of the `cmd` function which also runs
        several checks to alleviate debuggging TDEP issues. This
        cmd should not have `mpirun -np 4` or any other launcher.
        """
        pass

    @classmethod
    def input_files_present(cls, dir : PathLike) -> bool:
        dir = Path(dir) if not isinstance(dir, Path) else dir

        if not dir.is_dir():
            raise ValueError(f"The provided path '{dir}' is not a valid directory.")

        required_files = cls.input_files()
        return all((dir / file).exists() for file in required_files)
    
    def cmd(self, ncores, launcher = "mpirun"):
        """
        Checks if all reuqired input files are present and
        that the parameters passed to the cmd are valid. If so
        the command to run the program is generated and returned.
        """
        if self.input_files_present and self.input_parameters_valid():
            return self.__cmd()
        
    def mpirun(self, ncores : int, dir : str = "") -> int:
        initial_dir = os.getcwd()
        os.chdir(dir)
        res = os.system(f"mpirun -np {ncores} " + self.cmd())
        os.chdir(initial_dir)
        return res


