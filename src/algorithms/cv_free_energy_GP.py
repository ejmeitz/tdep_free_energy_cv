from typing import List
import matplotlib.pyplot as plt
import logging
import os, shutil
import numpy as np
from os.path import join
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel

from src import (
    # HeatCapFreeEnergyParams,
    AnharmonicFreeEnergy,
    write_tdep_meta,
    get_n_atoms_from_dump,
    remove_dump_headers,
    LammpsSimulator,
    temp_to_str
)

from .configs import HeatCapFreeEnergyParams, Paths

def init_F_arrs(N):
    U0 = np.zeros(N)
    F_ph = np.zeros(N)
    F_3 = np.zeros(N)
    F_4 = np.zeros(N)
    return U0, F_ph, F_3, F_4

def init_cv_arrs(N):
    return np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)

def plot_GP(GP, x_data, y_data, outpath):
    x_new = np.linspace(np.amin(x_data), np.amax(x_data), 1000);
    plt.plot(x_new, GP(x_new.reshape(-1, 1)), label='GP');
    plt.scatter(x_data, y_data, label = "Simulated Data");
    plt.legend();
    plt.xlabel("Tempearture");
    plt.ylabel("Irred IFC");
    plt.savefig(outpath);
    plt.close();

def fit_GP(x_data, y_data, plot_outpath):
    k = ConstantKernel(1.0) * RBF(length_scale=0.25, length_scale_bounds=(1e-2, 1e3)) +\
            WhiteKernel(1e-2, noise_level_bounds=(1e-7,1e4))
    gpr = GaussianProcessRegressor(kernel = k, n_restarts_optimizer=9, normalize_y=True).fit(x_data.reshape(-1, 1), y_data)
    return gpr

def get_fd_coeffs(size):
    if size == 3:
        return [1, -2, 1], [-1, 0, 1]
    elif size == 5:
        return [-1/12, 4/3, -5/2, 4/3, -1/12], [-2, -1, 0, 1, 2]
    elif size == 7:
        return [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90], [-3,-2,-1,0,1,2,3]
    else:
        raise ValueError(f"Only support stencils of size 3,5,7 for second-order central diff. Got {size}")

def cv_norm_from_F(F, coeffs, T, dT, kB = 8.61733262e-5):
    return -T * np.sum(F*coeffs) / (3*dT*dT*kB)

def plot_cv(Ts, cvs, outpath, title = None):
    plt.figure(dpi = 300);
    plt.hlines(1.0, np.amin(Ts), np.amax(Ts), color = "black", linestyles="dashed", lw = 1.5)
    plt.scatter(Ts, cvs, color = "#f76df7", s = 50);
    plt.xlabel("Temperature [K]");
    plt.ylabel("Cv / 3*N*$k_B$");
    if title is not None:
        plt.title(title);
    plt.savefig(outpath);
    

def run_cv_free_energy_GP(p : HeatCapFreeEnergyParams, paths : Paths):

    # Run free energy at all temperatures
    # Fit GP to the free energy and differentiate that
    free_energy_dir = join(paths.basepath, "FREE_ENERGY_CALCS")
    os.mkdir(free_energy_dir)
    results_dir = join(paths.basepath, "RESULTS")
    os.mkdir(results_dir)

    U0, F_ph, F3, F4 = init_F_arrs(len(p.temperatures))
    afe = AnharmonicFreeEnergy(p.k_mesh, quantum = p.quantum, stochastic = p.stochastic)
    for i,T in enumerate(p.temperatures):
        temp_str = temp_to_str(T)
        free_energy_dir_T = join(free_energy_dir, f"T{temp_str}")
        os.mkdir(free_energy_dir_T)

        shutil.copyfile(paths.ucposcar_path, join(free_energy_dir_T, "infile.ucposcar"))
        shutil.copyfile(paths.ssposcar_path, join(free_energy_dir_T, "infile.ssposcar"))
    
        # Re-use MD data that was used get the IFCs
        sim_dir = join(paths.ifc_path, "simulation_data", f"T{temp_str}")
        shutil.copyfile(join(sim_dir, "infile.positions"), join(free_energy_dir_T, "infile.positions"))
        shutil.copyfile(join(sim_dir, "infile.forces"), join(free_energy_dir_T, "infile.forces"))
        shutil.copyfile(join(sim_dir, "infile.meta"), join(free_energy_dir_T, "infile.meta"))
        shutil.copyfile(join(sim_dir, "infile.stat"), join(free_energy_dir_T, "infile.stat"))
        # shutil.copyfile(join(sim_dir, "outfile.U0"), join(free_energy_dir_T, "outfile.U0"))

        # Move necessary force constants
        shutil.copyfile(join(sim_dir, "outfile.forceconstant"), join(free_energy_dir_T, "infile.forceconstant"))
        shutil.copyfile(join(sim_dir, "outfile.forceconstant_thirdorder"),
                         join(free_energy_dir_T, "infile.forceconstant_thirdorder"))
        shutil.copyfile(join(sim_dir, "outfile.forceconstant_fourthorder"), 
                          join(free_energy_dir_T, "infile.forceconstant_fourthorder"))
        
        # Run free energy calculation
        afe.mpirun(p.n_cores, free_energy_dir_T)

        F_ph[i], F3[i], F4[i], _, _ = afe.parse_bulk_F_from_log(free_energy_dir_T)

        # Get U0 from prior MD simulation, its also in the log file but harder to parse
        U0[i] = np.loadtxt(join(sim_dir, "outfile.U0"))[-1]

    # Fit GP to each component and total
    F_total = U0 + F_ph + F3 + F4 # missing cumulants technically

    GP_total = fit_GP(p.temperatures, F_total)
    GP_U0 = fit_GP(p.temperatures, U0)
    GP_ph = fit_GP(p.temperatures, F_ph)
    GP_F3 = fit_GP(p.temperatures, F3)
    GP_F4 = fit_GP(p.temperatures, F4)

    # Plot GPs
    plot_GP(GP_total, p.temperatures, F_total, join(paths.basepath, "F_total_GP.png"))
    plot_GP(GP_U0, p.temperatures, U0, join(paths.basepath, "U0_GP.png"))
    plot_GP(GP_F3, p.temperatures, F3, join(paths.basepath, "F3_GP.png"))
    plot_GP(GP_F4, p.temperatures, F4, join(paths.basepath, "F4_GP.png"))
    plot_GP(GP_ph, p.temperatures, F_ph, join(paths.basepath, "F_ph_GP.png"))

    # Figure out temps needed for FD stencil
    coeffs, temp_offsets = get_fd_coeffs(p.fd_stencil_size)
    temp_stencil = np.array([to*p.dT for to in temp_offsets])

    cv_total, cv_ph, cv_F3, cv_F4, cv_U0 = init_cv_arrs(len(p.temperatures))
    for i,T in enumerate(p.temperatures):
        new_temps = T + temp_stencil
        # Evaluate GPs at each temperature on stencil
        F_total_stencil = GP_total.predict(new_temps.reshape(-1, 1))
        U0_stencil = GP_U0.predict(new_temps.reshape(-1, 1))
        F_ph_stencil = GP_ph.predict(new_temps.reshape(-1, 1))
        F3_stencil = GP_F3.predict(new_temps.reshape(-1, 1))
        F4_stencil = GP_F4.predict(new_temps.reshape(-1, 1))

        cv_total[i] = cv_norm_from_F(F_total_stencil, coeffs, T, p.dT)
        cv_U0[i] = cv_norm_from_F(U0_stencil, coeffs, T, p.dT)
        cv_F3[i] = cv_norm_from_F(F3_stencil, coeffs, T, p.dT)
        cv_F4[i] = cv_norm_from_F(F4_stencil, coeffs, T, p.dT)
        cv_ph[i] = cv_norm_from_F(F_ph_stencil, coeffs, T, p.dT)

        np.savetxt(join(results_dir, f"free_energies_{temp_str}.txt"), np.column_stack([new_temps, U0_stencil, F_ph_stencil, F3_stencil, F4_stencil, F_total_stencil]),
                    header = "# T U0 F_ph F_3 F_4 F_total", fmt = "%.12f")
    
    cv_sum = cv_U0 + cv_ph + cv_F3 + cv_F4
    np.savetxt(join(results_dir, f"heat_capacity_{temp_str}.txt"), np.column_stack([T, cv_U0, cv_ph, cv_F3, cv_F4, cv_sum, cv_total]),
                    fmt = "%.10f", header = "Cv / 3 N kB decomposition\nTemperature U0_part F_ph_part F3_part F4_part Cv_sum Cv_total")
        
    # Summary plots for heat capacity
    plot_cv(p.temperatures, cv_total, join(results_dir, "cv_total.png"))



        # Each free energy run should produce which we can get mode heat capacities from:
        # - outfile.F_ph_mode_resolved
        # - outfile.delta_F3_mode_resolved
        # - outfile.delta_F4_mode_resolved

        # Check convergence of U0 and cum2 somehow??



