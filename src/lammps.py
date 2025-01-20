import os, shutil, subprocess
from typing import List, Optional
import numpy as np

from .util import PathLike

def get_n_atoms_from_dump(dump_path : PathLike):
    try:
        # Run the bash command
        result = subprocess.run(
            f"head -n 4 {dump_path} | tail -n 1",
            shell=True,
            check=True,
            text=True,
            stdout=subprocess.PIPE
        )
        # Extract the output and convert to an integer
        output_line = result.stdout.strip()
        return int(output_line)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
    except ValueError as e:
        print(f"Could not convert output to integer: {e}")

def remove_dump_headers(simulation_folder):
        bash_script = f""" #!/bin/bash

            # figure out how many atoms there are
            na=`head -n 4 dump.forces | tail -n 1`
            # remove the header from the stat file
            grep -v '^#' dump.stat > ./infile.stat
            # figure out how many timesteps there are
            nt=`wc -l ./infile.stat"""
        
        bash_script_end = """| awk '{print $1}'`
        
            # create the positions and force files
            [ -f infile.forces ] && rm infile.forces
            [ -f infile.positions ] && rm infile.positions
            for t in `seq 1 ${nt}`
            do
                nl=$(( ${na}+9))
                nll=$(( ${nl}*${t} ))
                head -n ${nll} dump.forces | tail -n ${na} >> infile.forces
                head -n ${nll} dump.positions | tail -n ${na} >> infile.positions
            done
        """

        bash_script += bash_script_end

        bash_file = os.path.join(simulation_folder, "remove_dump_headers")
        with open(bash_file, "w") as f:
            f.write(bash_script)
        
        os.chdir(simulation_folder)
        os.system(f"chmod u+x {bash_file}")
        os.system(f"{bash_file}")



class InFile():

    '''
    Class to keep track of, and modify variables for an in-file
    '''

    def __init__(self, path, required_keys : Optional[List[str]] = None):
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

        if required_keys is not None:
            keys_present = np.array([k in self.free_variables for k in required_keys])
            if not all(keys_present):
                raise ValueError(f"Missing required variables in LAMMPS file: {required_keys[~keys_present]}")
        

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
                        else:
                            raise RuntimeError(f"Error parsing LAMMPS infile, got unknown mode : {mode}")
                        
                        self.free_variable_modes[var] = mode
                        self.free_variable_line_numbers[var] = line_number
                    if mode == "string":
                        self.free_variables[var] = value
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

# All of these could easily be combined into one that just takes a dict 
# with the keys and values as a param.
class LammpsSimulator:
    def __init__(
            self,
            base_infile_path : PathLike,
            run_dir : PathLike,
            edited_variables_dict : dict,
            required_keys : Optional[List[str]] = None
        ):
        self.base_infile_path = base_infile_path
        self.run_dir = run_dir

        # Copy infile to run dir
        base_infile_name = os.path.basename(self.base_infile_path)
        self.new_infile_path = os.path.join(run_dir, base_infile_name)
        shutil.copyfile(self.base_infile_path, self.new_infile_path)

        # Modify variables of in-file
        required_keys.append("velocity_seed")
        self.infile = InFile(self.new_infile_path, required_keys)
        self.infile.edit_variables(edited_variables_dict)
        self.infile.edit_variables({"velocity_seed" : np.random.randint(1000,1000000)})

    def run(self, dir : PathLike) -> int:
        os.chdir(dir)
        return os.system(f"lmp -in {self.new_infile_path} -screen lmp.screen")
    
    def mpirun(self, dir : PathLike, ncores : int) -> int:
        os.chdir(dir)
        return os.system(f"mpirun -np {ncores} lmp -in {self.new_infile_path} -screen lmp.screen")
