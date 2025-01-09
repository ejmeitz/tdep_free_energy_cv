import logging
import os
from pathlib import Path
from typing import Union

PathLike = Union[Path, str]

def setup_logging(logname : str, savepath : PathLike = os.getcwd()):
    logging.basicConfig(
        level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(savepath, logname)), 
            logging.StreamHandler()
        ]
    )
    logging.info("Logger initialized")


import os

class LAMMPSInFile():

    '''
    Class to keep track of, and modify variables for an in-file
    '''

    def __init__(self, path):
        '''
        Path: Absolute file path to the in-file, contains filename and extension.
        e.g. C:/Users/ejmei/Desktop/in.argon
        '''
        self.path = path

        '''
        Name: Filename with extension
        Free Variables: All variables which are purely numeric. Any variable that is a function
            of other variables or a string is ignored.
        Free Variable Line Number: 0 indexed line number corresponding to the line in the in-file
            where that variable is found
        '''
        self.name = os.path.basename(path)

        self.free_variables = {}
        self.free_variable_line_numbers = {}
        self.free_variable_modes = {}
        self.__parse_variables()
        

    def __parse_variables(self) -> None:
        '''
        Opens in-file and extracts parameters which are specified as 'variables' AND 
            do not depend on other variables or output from LAMMPS.

        Example: 
            "variable dt equal 0.002" will be parsed
            "variable t_damp equal 100*dt" will NOT be parsed
        '''

        with open(self.path, 'r') as f:
            for line_number, line in enumerate(f.readlines()):
                if line.strip().startswith("variable"):
                    _, var, mode, value, *_ = line.strip().split()
                    if value.replace('.', '', 1).isdigit():
                        if mode == "equal":
                            self.free_variables[var] = float(value)
                        elif mode == "string":
                            self.free_variables[var] = value
                        else:
                            raise RuntimeError(f"Error parsing LAMMPS infile, got unknown mode : {mode}")
                        
                        self.free_variable_modes[var] = mode
                        self.free_variable_line_numbers[var] = line_number

    
    def edit_variables(self, changed_variables : dict) -> None:
        '''
        Opens the in-file for this job and modifies the variable values to match those
        pass in the constructor.

        Changed Variables: A dictionary where the keys are the variable to change and the
            value is the updated value to write to the in-file
        '''
        #Copy file contents
        with open(self.path, 'r+') as in_file:
            lines = in_file.readlines()
        #Edit file contents
        for new_variable, value in changed_variables.items():
            if new_variable in self.free_variables.keys():
                idx = self.free_variable_line_numbers[new_variable]
                mode = self.free_variable_modes[new_variable]
                lines[idx] = f"variable {new_variable} {mode} {value}\n"
            else:
                print("============================================================================")
                print(f"{new_variable} is not a modifiable variable in the in-file. This variable will not be changed. In-file located at: {self.path}")
                print("============================================================================")
        #Re-write file contents
        with open(self.path,'w') as in_file:
            in_file.writelines(lines)