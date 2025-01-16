from dataclasses import dataclass
from typing import Optional, List
import logging
import os, shutil
import numpy as np
from os.path import join

from src import (
    PathLike,
    InterpolateIFCParams,
    AnharmonicFreeEnergy,
    write_tdep_meta,
    get_n_atoms_from_dump,
    run_interpolate_irred,
    remove_dump_headers,
    LammpsSimulator
)

@dataclass
class HeatCapFreeEnergyParams:
    temperature : float
    dT : float
    k_mesh : List[int]
    ucposrcar_path : PathLike
    ssposcar_path : PathLike
    interp_settings : InterpolateIFCParams
    n_cores : int = 1
    fd_stencil_size : int = 3
    quantum : bool = False # use Bose-Einstein or Classical Occupation
    stochastic : bool = False # if 2nd order IFCs are fit with sTDEP then True
    force_calc : str = "lammps"
    script_name : Optional[str] = None #(e.g., LJ_argon_dynamics.in)
    structure_path : Optional[str] = None #(e.g. /home/emeitz/initial_structure_LJ.data)

# Heat capacity is second derivative of free energy
def get_fd_coeffs(size):
    if size == 3:
        return [1, -2, 1], [-1, 0, 1]
    elif size == 5:
        return [-1/12, 4/3, -5/2, 4/3, -1/12], [-2, -1, 0, 1, 2]
    elif size == 7:
        return [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90], [-3,-2,-1,0,1,2,3]
    else:
        raise ValueError(f"Only support stencils of size 3,5,7 for second-order central diff. Got {size}")

def check_params_consistent(p : HeatCapFreeEnergyParams, temp_stencil : List[float]) -> None:

    if not set(temp_stencil).issubset(set(p.interp_settings.temps_to_interpolate)):
        raise ValueError("Interpolation settings inconsistent with FD stencil. Not all temps on stencil calculated by interpolation.")

    if p.interp_settings.rc3 is None or p.interp_settings.rc4 is None:
        raise ValueError("Missing anharmonic IFC cutoff, heat capacity calculation requires 2nd, 3rd and 4th order IFCs.")

    if p.interp_settings.make_ss_ifcs == False:
        logging.warning("Setting `make_ss_ifcs` to True in interp_settings. This is required for heat capacity calculation.")
        p.interp_settings.make_ss_ifcs = True

    if p.interp_settings.interpolate_U0 == False:
        logging.wraning("Setting `interpolate_U0` to True in interp settings. This is required for heat capacity calculation.")
        p.interp_settings.interpolate_U0 = True

    if p.n_cores == 1:
        logging.warning("Got 1 core. You should really parallelize the anharmonic free energy calculation...")

def run_dummy_lammps(T, structure_path, base_infile_path, sim_root_dir):
    N_steps = 10_000 # just something small this data isnt used
    data_interval = 5_000
    N_configs = 3 #(0, 5000, 10000)

    var_dict = {"T" : T, "N_steps" : N_steps, "data_interval" : data_interval, "structure_path" : structure_path}
    required_vars = ["T", "N_steps", "structure_path", "data_interval"]
    ls = LammpsSimulator(base_infile_path, sim_root_dir, var_dict, required_vars)

    ls.run(sim_root_dir)

    N_atoms = get_n_atoms_from_dump(join(sim_root_dir, "dump.positions"))
    remove_dump_headers(sim_root_dir) # makes infile.positions, infile.stat, infile.forces
    write_tdep_meta(sim_root_dir, N_atoms, N_configs, 1.0, T) #timestep arbitrary

    return N_configs, N_atoms

def run_cv_free_energy(p : HeatCapFreeEnergyParams):


    coeffs, temp_offsets = get_fd_coeffs(p.fd_stencil_size)
    temp_stencil = p.temperature + np.array([to*p.dT for to in temp_offsets])

    check_params_consistent(p, temp_stencil)

    # Run interpolation, expects infile.ucposcar and infile.ssposcar in root dir
    shutil.copyfile(p.ucposrcar_path, join(p.lds.basepath, "infile.ucposcar"))
    shutil.copyfile(p.ssposcar_path, join(p.lds.basepath, "infile.ssposcar"))
    run_interpolate_irred(p.interp_settings)

    # Generate required input files at central temperature
    # these are just dummy files to we can call anharmonic_free_energy
    md_dir = join(p.basepath, "CENTRAL_TEMP_SIMULATION")
    os.mkdir(md_dir)
    if p.force_calc == "lammps":
        logging.info(f"Running LAMMPS at {T} K")

        file_path = os.path.dirname(os.path.realpath(__file__))
        base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.script_name)
        base_infile_path = os.path.abspath(base_infile_path)

        N_configs, N_atoms = run_dummy_lammps(p.structure_path, T, base_infile_path, md_dir)

    # Run free energy calculation on each of the interpolation points
    free_energy_dir = join(p.basepath, "FREE_ENERGY_CALCS")
    os.mkdir(free_energy_dir)

    afe = AnharmonicFreeEnergy(p.k_mesh, quantum = p.quantum, stochastic = p.stochastic)
    for T in temp_stencil:
        free_energy_dir_T = join(free_energy_dir, f"T{T}")
        os.mkdir(free_energy_dir_T)

        shutil.copyfile(p.ucposrcar_path, join(free_energy_dir_T, "infile.ucposcar"))
        shutil.copyfile(p.ssposcar_path, join(free_energy_dir_T, "infile.ssposcar"))

        # Move MD simulation data to this dir
        # This will technically be at the wrong temperature.
        # This only affects the calculation of U0, so we will 
        # ignore U0 calculated by this and use the interpolated value.
        shutil.copyfile(join(md_dir, "infile.stat"), join(free_energy_dir_T, "infile.stat"))
        shutil.copyfile(join(md_dir, "infile.forces"), join(free_energy_dir_T, "infile.forces"))
        shutil.copyfile(join(md_dir, "infile.positions"), join(free_energy_dir_T, "infile.positions"))

        # The infile.meta MUST have the correct temperature
        # as it is used to get phonon occupations occupations
        write_tdep_meta(free_energy_dir_T, N_atoms, N_configs, 1.0, T)

        afe.mpirun(p.n_cores, free_energy_dir_T)


    # Each free energy run should produce which we can get mode heat capacities from:
    # - outfile.F_ph_mode_resolved
    # - outfile.delta_F3_mode_resolved
    # - outfile.delta_F4_mode_resolved
    # Use the interpolated U0s to account for its affect on heat capacity
    # Parse the bulk free energy from default TDEP and compare with my results



