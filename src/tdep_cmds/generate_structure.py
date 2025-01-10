from typing import List

from .tdep_cmd import TDEP_Command

class GenerateStructure(TDEP_Command):

    def __init__(
            self,
            dims : List[int],
            log_file : str = "generate_structure.log"
        ):
        self.dims = dims
        self.log_file = log_file

    @property
    def input_files(self):
        return ["infile.ucposcar"]
    
    def output_files(self):
        return ["outfile.ssposcar"]
    
    def input_parameters_valid(self):
        if len(self.dims) != 3:
            raise RuntimeError("Super cell dims must have length 3")
        return True

    def _cmd(self) -> str:
        return f"generate_structure -d {self.dims[0]} {self.dims[1]} {self.dims[2]} > {self.log_file}"
    