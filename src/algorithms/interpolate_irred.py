import numpy as np
from numpy.polynomial.polynomial import Polynomial
import scipy.interpolate as interp
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
import logging
import os, shutil
from os.path import join
import matplotlib.pyplot as plt

from src import (
    ExtractForceConstants,
    write_tdep_meta,
    temp_to_str
)

from .ifc_from_MD import run_ifc_from_MD
from .configs import InterpolateIFCParams, IFC_MD_Params, Paths


def interpolate(mode, X_INTERP, X_DATA, Y_DATA):

    if mode == "linear":
        return np.interp(X_INTERP, X_DATA, Y_DATA)
    elif mode == "cubic":
        cs = interp.CubicSpline(X_DATA, Y_DATA)
        return cs, cs(X_INTERP)
    elif mode == "gpr": # not really interpolation but whatever
        k = ConstantKernel(1.0) * RBF(length_scale=0.1, length_scale_bounds=(1e-2, 1e4)) +\
              WhiteKernel(1e-2, noise_level_bounds=(1e-7,1e3))
        gpr = GaussianProcessRegressor(kernel = k, n_restarts_optimizer=9, normalize_y=True).fit(X_DATA.reshape(-1, 1), Y_DATA)
        return gpr.predict, gpr.predict(X_INTERP.reshape(-1, 1))
        
    # elif mode == "lagrange": #numerically unstable
    #     poly = Polynomial(interp.lagrange(X_DATA, Y_DATA).coef[::-1])
    #     return poly, poly(X_INTERP)
    else:
        raise ValueError(f"Invalid interpolate mode : {mode}. Expected linear or cubic.")

def simulate_interp_node_data(p : InterpolateIFCParams, paths : Paths):
    # Get force constants
    if p.force_calc == "lammps":
        # We will clean up manually at end, might need this sim data later
        ifc_params = IFC_MD_Params(p.temps_to_simulate, p.n_cores_max, p.rc2, p.rc3, p.rc4, p.lds, False)
        sim_dir, ifc_dir, mean_sim_temps = run_ifc_from_MD(ifc_params, paths)
        logging.info("Using mean simulation temperatures for interpolation.")
    elif p.force_calc == "vasp":
        raise NotImplementedError("VASP force calculator not implemented yet")
    else:
        raise ValueError(f"Unknown force calculator, expected lammps or vasp, got : {p.force_calc}")
    
    return sim_dir, ifc_dir, mean_sim_temps

def plot_interpolant(x_sim, y_sim, interp_poly : Polynomial, outpath):
    x_new = np.linspace(np.amin(x_sim), np.amax(x_sim), 1000);
    plt.plot(x_new, interp_poly(x_new.reshape(-1, 1)), label='Interpolant');
    plt.scatter(x_sim, y_sim, label = "Simulated Data");
    plt.legend();
    plt.xlabel("Tempearture");
    plt.ylabel("Irred IFC");
    plt.savefig(outpath);
    plt.close();


def run_interpolate_irred(p : InterpolateIFCParams, paths : Paths,
                             sim_dir, ifc_dir, mean_sim_temps):

    # Check that interp temps are not actually extrapolating
    min_T_sim = np.amin(p.simulated_temps)
    max_T_sim = np.amax(p.simulated_temps)
    min_T_interp = np.amin(p.temps_to_interpolate)
    max_T_interp = np.amax(p.temps_to_interpolate)

    if max_T_interp > max_T_sim or min_T_interp < min_T_sim:
        raise ValueError("Temperatures passed require extrapolation, exiting.")
    
    TEMPS_SIM = np.sort(p.simulated_temps)
    TEMPS_INTERP = np.sort(p.temps_to_interpolate)


    # Parse irreducible force constants
    irred2_ifcs = []; irred3_ifcs = []; irred4_ifcs = []; U0s = []
    for i, T in enumerate(TEMPS_SIM):
        temp_str = temp_to_str(T)
        ifc_dir_T = join(ifc_dir, f"T{temp_str}")

        # Store irred ifcs
        irred2_ifcs.append(np.loadtxt(join(ifc_dir_T, f"infile.irrifc_secondorder")))

        if p.rc3 is not None:
            irred3_ifcs.append(np.loadtxt(join(ifc_dir_T, f"infile.irrifc_thirdorder")))
        
        if p.rc4 is not None:
            irred4_ifcs.append(np.loadtxt(join(ifc_dir_T, f"infile.irrifc_fourthorder")))

        if p.interpolate_U0:
            tmp = np.loadtxt(join(ifc_dir_T, "outfile.U0"))
            U0s.append(tmp[-1])

    # Interpolate the irreducible IFCs
    irred2_ifcs = np.array(irred2_ifcs)
    new_irred2_ifcs = np.zeros((irred2_ifcs.shape[1], len(TEMPS_INTERP)))
    for i in range(irred2_ifcs.shape[1]):
        func, new_irred2_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred2_ifcs[:,i])

        #* TODO JUST PLOT THEM ALL....
        if i == 0:
            plot_interpolant(TEMPS_SIM, irred2_ifcs[:,0], func, join(paths.basepath, "irred_interp_SO.png"))

    if p.rc3 is not None:
        irred3_ifcs = np.array(irred3_ifcs)
        new_irred3_ifcs = np.zeros((irred3_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred3_ifcs.shape[1]):
            func, new_irred3_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred3_ifcs[:,i])

            if i == 0:
                plot_interpolant(TEMPS_SIM, irred3_ifcs[:,0], func, join(paths.basepath, "irred_interp_TO.png"))

    if p.rc4 is not None:
        irred4_ifcs = np.array(irred4_ifcs)
        new_irred4_ifcs = np.zeros((irred4_ifcs.shape[1], len(TEMPS_INTERP)))
        for i in range(irred4_ifcs.shape[1]):
            func, new_irred4_ifcs[i,:] = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, irred4_ifcs[:,i])

            if i == 0:
                plot_interpolant(TEMPS_SIM, irred4_ifcs[:,0], func, join(paths.basepath, "irred_interp_FO.png"))

    if p.interpolate_U0:
        func, interpolated_U0s = interpolate(p.interp_mode, TEMPS_INTERP, mean_sim_temps, U0s)
        plot_interpolant(TEMPS_SIM, U0s, func, join(paths.basepath, "irred_interp_U0.png"))
        


    # Save interpolated IFCs per tempearture
    irred2_paths = []
    irred3_paths = [] if p.rc3 is not None else None
    irred4_paths = [] if p.rc4 is not None else None
    irred_out_path = join(paths.basepath, "INTERPOLATED_IRRED_IFCS")
    os.mkdir(irred_out_path)
    for j, T in enumerate(TEMPS_INTERP):
        temp_str = temp_to_str(T)
        p2 = join(irred_out_path, f"infile.irrifc_secondorder_{temp_str}")
        np.savetxt(p2, new_irred2_ifcs[:,j])
        irred2_paths.append(p2)

        if p.rc3 is not None:
            p3 = join(irred_out_path, f"infile.irrifc_thirdorder_{temp_str}")
            np.savetxt(p3, new_irred3_ifcs[:,j])
            irred3_paths.append(p3)

        if p.rc4 is not None:
            p4 = join(irred_out_path, f"infile.irrifc_fourthorder_{temp_str}")
            np.savetxt(p4, new_irred4_ifcs[:,j])
            irred4_paths.append(p4)
    
    if p.interpolate_U0:
        np.savetxt(join(irred_out_path, "interpolated_U0s.txt"), np.column_stack([TEMPS_INTERP,interpolated_U0s]),
                     fmt = "%.12f", header = "# Temperature U0")
    #* TODO INTERPOLATE SECOND ORDER CUMULANT

    # Writes out command that must be used to get the same forcemap
    ef_tmp = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    with open(join(paths.basepath, "READ_IRRED_IFC_CMD.txt"), "w") as f:
        f.write("""# YOU MUST USE THIS COMMAND WHEN USING THE INTERPOLATED IFCS.
                # THE PARAMS MUST MATCH THOSE USED TO CREATE THE INTERPOLATION POINTS.
                # THE infile.ucposcar AND infile.ssposcar MUST ALSO BE IDENTICAL."""
                )
        f.write(ef_tmp._cmd())


    if p.make_ss_ifcs:
        logging.info("Converting interpolated irreducible IFCs, back to super cell IFCs")
        ss_out_dir = convert_irred_to_ss_ifc(p, paths, TEMPS_INTERP, sim_dir, irred_out_path)
        shutil.copyfile(join(irred_out_path, "interpolated_U0s.txt"), join(ss_out_dir, "interpolated_U0s.txt"))

    if p.cleanup:
        os.chdir(paths.basepath)
        os.system(f"rm -rf {sim_dir}")

    return irred_out_path, ss_out_dir
    
def convert_irred_to_ss_ifc(p : InterpolateIFCParams, paths : Paths,
                             TEMPS_INTERP, sim_dir, irred_out_path):
    
    ss_out_dir = join(paths.basepath, "INTERPOLATED_SUPERCELL_IFCS")
    os.mkdir(ss_out_dir)

    ef = ExtractForceConstants(p.rc2, p.rc3, p.rc4, read_irreducible=True)
    for T in TEMPS_INTERP:
        temp_str = temp_to_str(T)
        sim_dir_T = join(sim_dir, f"T{temp_str}_SS_RECONSTRUCT")
        os.mkdir(sim_dir_T)
        os.chdir(sim_dir_T)

        # Copy interpolated irreducible IFCs to this dir
        shutil.copyfile(join(irred_out_path, f"infile.irrifc_secondorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_secondorder"))
        if p.rc3 is not None:
            shutil.copyfile(join(irred_out_path, f"infile.irrifc_thirdorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_thirdorder"))
        if p.rc4 is not None:
            shutil.copyfile(join(irred_out_path, f"infile.irrifc_fourthorder_{temp_str}"), join(sim_dir_T, "infile.irrifc_fourthorder"))

        # Move required files
        shutil.copyfile(join(paths.basepath, "infile.ucposcar"), join(sim_dir_T, "infile.ucposcar"))
        shutil.copyfile(join(paths.basepath, "infile.ssposcar"), join(sim_dir_T, "infile.ssposcar"))

        # forces, positions and stat file will not affect
        # IFC calculation but need to be valid.
        # Make some fake ones
        # shutil.copyfile(join(sim_dir, f"T{TEMPS_SIM[0]}", "infile.positions"), join(sim_dir_T, "infile.positions"))
        # shutil.copyfile(join(sim_dir, f"T{TEMPS_SIM[0]}", "infile.forces"), join(sim_dir_T, "infile.forces"))
        ssposcar_positions = np.loadtxt(join(paths.basepath, "infile.ssposcar"), skiprows = 8)
        N_atoms = ssposcar_positions.shape[0]
        np.savetxt(join(sim_dir_T, "infile.positions"), ssposcar_positions, fmt = "%.7f")
        np.savetxt(join(sim_dir_T, "infile.forces"), np.zeros((N_atoms, 3)), fmt = "%.7f")
        stat_fmt = ['%d', '%d'] + ['%.8f'] * 11
        np.savetxt(join(sim_dir_T, "infile.stat"), np.zeros((1, 13)), fmt = stat_fmt) # this isnt even used, but still has to exist 
        write_tdep_meta(sim_dir_T, N_atoms, 1, 1.0, T)

        # Re-run extract force constants using the irreducible IFCs
        # This still needs infile.meta, infile.stat, infile.positions and infile.forces
        # to run although those will not be used for calculating the IFCs
        ef.mpirun(p.n_cores_max, sim_dir_T)

        # Move force constants to make folder to avoid clean up 
        shutil.copyfile(join(sim_dir_T, "outfile.forceconstant"), join(ss_out_dir, f"infile.forceconstant_{temp_str}"))
        
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_thirdorder"), join(ss_out_dir, f"infile.forceconstant_thirdorder_{temp_str}"))
            
        if p.rc3 is not None:
            shutil.copyfile(join(sim_dir_T, "outfile.forceconstant_fourthorder"), join(ss_out_dir, f"infile.forceconstant_fourthorder_{temp_str}"))

    return ss_out_dir
