import os, shutil

from .util import PathLike


class InFile():

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

class LammpsCalculator:

    def __init__(
            self,
            base_infile_path : PathLike,
            run_dir : PathLike,
            structure_pattern : str, # expects filename w/ numbers removed (e.g. contcar_conf)
            n_configs : int,
            *,
            structure_key : str = "structure_path",
            config_key : str = "n_configs"
        ):
        self.base_infile_path = base_infile_path
        self.structure_pattern = structure_pattern
        self.run_dir = run_dir

        # Copy infile to run dir
        base_infile_name = self.base_infile_path.basename()
        self.new_infile_path = os.path.join(run_dir, base_infile_name)
        shutil.copyfile(self.base_infile_path, self.new_infile_path)

        # Modify variables of in-file
        self.infile = InFile(self.new_infile_path)
        self.infile.edit_variables({structure_key : os.path.join(run_dir, structure_pattern),
                                    config_key : n_configs})

    def run(self) -> int:
        return os.system(f"lmp < {self.new_infile_path}")
    
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
                echo "t ${t} ${nl} ${nll}"
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



