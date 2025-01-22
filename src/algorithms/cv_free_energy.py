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
    LammpsSimulator,
    temp_to_str
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

def check_params_consistent(p : HeatCapFreeEnergyParams) -> None:

    if p.n_cores == 1:
        logging.warning("Got 1 core. You should really parallelize the anharmonic free energy calculation...")

    if p.interp_settings.rc3 is None or p.interp_settings.rc4 is None:
        raise ValueError("Missing anharmonic IFC cutoff, heat capacity calculation requires 2nd, 3rd and 4th order IFCs.")

    if p.interp_settings.make_ss_ifcs == False:
        logging.warning("Setting `make_ss_ifcs` to True in interp_settings. This is required for heat capacity calculation.")
        p.interp_settings.make_ss_ifcs = True

    if p.interp_settings.interpolate_U0 == False:
        logging.wraning("Setting `interpolate_U0` to True in interp settings. This is required for heat capacity calculation.")
        p.interp_settings.interpolate_U0 = True

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

def run_cv_free_energy(p : HeatCapFreeEnergyParams, paths : Paths,
                        sim_dir, ifc_dir, mean_sim_temps):

    data_outpath = join(paths.basepath, "HeatCapData")
    os.mkdir(data_outpath)

    check_params_consistent(p)

    coeffs, temp_offsets = get_fd_coeffs(p.fd_stencil_size)
    INTERP_TEMPS = np.sort([ct + (to*p.dT) for to in temp_offsets for ct in p.temperatures])
    p.interp_settings.temps_to_interpolate = INTERP_TEMPS
    p.interp_settings.n_cores = p.n_cores

    # Run interpolation at all the necessary temps
    _, ss_out_dir = run_interpolate_irred(p.interp_settings, paths, sim_dir, ifc_dir, mean_sim_temps)


    for CENTER_TEMP in p.temperatures:
        temp_str = temp_to_str(CENTER_TEMP)
        data_outpath_T = join(data_outpath, temp_str)
        os.mkdir(data_outpath_T)

        # Generate required input files at central temperature
        # these are just dummy files so we can call anharmonic_free_energy
        # Needed to get U0, but we got that from interpolation
        md_dir = join(data_outpath_T, "DUMMY_SIMULATION")
        os.mkdir(md_dir)
        if p.force_calc == "lammps":
            logging.info(f"Running LAMMPS at {CENTER_TEMP} K")
            N_configs, N_atoms = run_dummy_lammps(CENTER_TEMP, p.interp_settings.lds.structure_path,
                                                p.interp_settings.lds.infile_path, md_dir)

        # Get interpolated U0 values
        U0_interpolated = np.loadtxt(join(ss_out_dir, "interpolated_U0s.txt")) # (temp, U0)
        U0_temps = U0_interpolated[:,0]; U0_values = U0_interpolated[:,1]
        U0_values = U0_values[np.argsort(U0_temps)] # just make sure in order

        # Run free energy calculation on each of the interpolation points
        free_energy_dir = join(data_outpath_T, "FREE_ENERGY_CALCS")
        os.mkdir(free_energy_dir)

        F_ph, F_3, F_4, cum2, cum3 = init_out_arrs(p.fd_stencil_size)
        afe = AnharmonicFreeEnergy(p.k_mesh, quantum = p.quantum, stochastic = p.stochastic)
        for i,T in enumerate(INTERP_TEMPS):
            temp_str = temp_to_str(T)
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
        #sgn = -1.0 if afe.stochastic else 1.0
        F_total = U0_values + F_ph + F_3 + F_4 #+ sgn*cum2 #+ cum3 # cum3 takes too long to converge to bother including

        for i in range(0,len(F_total) + 1, p.fd_stencil_size):
            s = slice(i,i+p.fd_stencil_size)
            bulk_cv = cv_norm_from_F(F_total[s], coeffs, CENTER_TEMP, p.dT)
            U0_deriv = cv_norm_from_F(U0_values[s], coeffs, CENTER_TEMP, p.dT)
            F_ph_deriv = cv_norm_from_F(F_ph[s], coeffs, CENTER_TEMP, p.dT)
            F3_deriv = cv_norm_from_F(F_3[s], coeffs, CENTER_TEMP, p.dT)
            F4_deriv = cv_norm_from_F(F_4[s], coeffs, CENTER_TEMP, p.dT)
            np.savetxt(join(data_outpath_T, "heat_capacity.txt"), np.column_stack([CENTER_TEMP, U0_deriv, F_ph_deriv, F3_deriv, F4_deriv, bulk_cv]),
                        fmt = "%.10f", header = "Cv / N kB decomposition\nTemperature U0_part F_ph_part F3_part F4_part Bulk_CV")

            np.savetxt(join(data_outpath_T, "free_energies.txt"), np.column_stack([INTERP_TEMPS[s], U0_values[s], F_ph[s], F_3[s], F_4[s], F_total[s]]),
                        header = "# T U0 F_ph F_3 F_4 F_total", fmt = "%.12f")

        # Each free energy run should produce which we can get mode heat capacities from:
        # - outfile.F_ph_mode_resolved
        # - outfile.delta_F3_mode_resolved
        # - outfile.delta_F4_mode_resolved

        # Check convergence of U0 and cum2 somehow??



