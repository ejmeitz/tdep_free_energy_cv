from typing import List
import logging
import os, shutil
import numpy as np
from os.path import join

from src import (
    # HeatCapFreeEnergyParams,
    AnharmonicFreeEnergy,
    write_tdep_meta,
    get_n_atoms_from_dump,
    remove_dump_headers,
    LammpsSimulator
)

from .interpolate_irred import run_interpolate_irred
from .configs import HeatCapFreeEnergyParams, Paths


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

def init_out_arrs(N):
    F_ph = np.zeros(N)
    F_3 = np.zeros(N)
    F_4 = np.zeros(N)
    cum2 = np.zeros(N)
    cum3 = np.zeros(N)
    return F_ph, F_3, F_4, cum2, cum3

# from TDEP F has units of eV/atom
def cv_norm_from_F(F, coeffs, T, dT, kB = 8.61733262e-5):
    return -T * np.sum(F*coeffs) / (3*dT*dT*kB)

def run_cv_free_energy(p : HeatCapFreeEnergyParams, paths : Paths):


    coeffs, temp_offsets = get_fd_coeffs(p.fd_stencil_size)
    temp_stencil = p.temperature + np.array([to*p.dT for to in temp_offsets])

    check_params_consistent(p, temp_stencil)
    p.interp_settings.n_cores = p.n_cores

    irred_out_path, ss_out_dir = run_interpolate_irred(p.interp_settings, paths)

    # Generate required input files at central temperature
    # these are just dummy files to we can call anharmonic_free_energy
    md_dir = join(paths.basepath, "DUMMY_SIMULATION")
    os.mkdir(md_dir)
    if p.force_calc == "lammps":
        logging.info(f"Running LAMMPS at {p.temperature} K")

        file_path = os.path.dirname(os.path.realpath(__file__))
        base_infile_path = os.path.join(file_path, "..", "..", "data", "lammps_scripts", p.interp_settings.lds.script_name)
        base_infile_path = os.path.abspath(base_infile_path)

        N_configs, N_atoms = run_dummy_lammps(p.temperature, p.interp_settings.lds.structure_path, base_infile_path, md_dir)

    # Get interpolated U0 values
    U0_interpolated = np.loadtxt(join(ss_out_dir, "interpolated_U0s.txt")) # (temp, U0)
    U0_temps = U0_interpolated[:,0]; U0_values = U0_interpolated[:,1]
    U0_values = U0_values[np.argsort(U0_temps)] # just make sure in order

    # Run free energy calculation on each of the interpolation points
    free_energy_dir = join(paths.basepath, "FREE_ENERGY_CALCS")
    os.mkdir(free_energy_dir)

    F_ph, F_3, F_4, cum2, cum3 = init_out_arrs(p.fd_stencil_size)
    afe = AnharmonicFreeEnergy(p.k_mesh, quantum = p.quantum, stochastic = p.stochastic)
    temps = np.sort(temp_stencil)
    for i,T in enumerate(temps):
        temp_str = f"{T}".replace('.', '_')
        free_energy_dir_T = join(free_energy_dir, f"T{temp_str}")
        os.mkdir(free_energy_dir_T)

        shutil.copyfile(paths.ucposcar_path, join(free_energy_dir_T, "infile.ucposcar"))
        shutil.copyfile(paths.ssposcar_path, join(free_energy_dir_T, "infile.ssposcar"))

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

        # Move force constants
        shutil.copyfile(join(ss_out_dir, f"infile.forceconstant_{temp_str}"),
                        join(free_energy_dir_T, "infile.forceconstant"))
        shutil.copyfile(join(ss_out_dir, f"infile.forceconstant_thirdorder_{temp_str}"),
                        join(free_energy_dir_T, "infile.forceconstant_thirdorder"))
        shutil.copyfile(join(ss_out_dir, f"infile.forceconstant_fourthorder_{temp_str}"),
                        join(free_energy_dir_T, "infile.forceconstant_fourthorder"))

        afe.mpirun(p.n_cores, free_energy_dir_T)

        F_ph[i], F_3[i], F_4[i], cum2[i], cum3[i] =\
            afe.parse_bulk_F_from_log(free_energy_dir_T)

    #* TODO CALCULATE 2nd ORDER CUMULANT IN INTERPOLATE IRRED
    # Cannot trust cumulant / U0 from the MD simulation
    #sgn = -1.0 if afe.stochastic else 1.0
    F_total = U0_values + F_ph + F_3 + F_4 #+ sgn*cum2 #+ cum3 # cum3 takes too long to converge to bother including

    cv = cv_norm_from_F(F_total, coeffs, p.temperature, p.dT)
    np.savetxt(join(paths.basepath, "heat_capacity.txt"), [cv], fmt = "%.12f")

    np.savetxt(join(paths.basepath, "free_energies.txt"), np.column_stack([temps, U0_values, F_ph, F_3, F_4, F_total]),
                header = "T U0 F_ph F_3 F_4 F_total", fmt = "%.12f")

    # Each free energy run should produce which we can get mode heat capacities from:
    # - outfile.F_ph_mode_resolved
    # - outfile.delta_F3_mode_resolved
    # - outfile.delta_F4_mode_resolved

    # Check convergence of U0 and cum2 somehow??



