import argparse
import logging
import os, shutil
import yaml

from src import (
    setup_logging,
    Paths,
    HeatCapFreeEnergyParams,
    run_cv_free_energy,
    initialize_free_energy, 
    run_interpolate_irred,
    InterpolateIFCParams,
    LammpsDynamicsSettings
)

"""
Calculates heat capacity from the second derivatve of free energy. Force constants are interpolated
to each temperature on the finite difference stencil. By default the interpolation is done at each call
but prior data can be used by setting `Paths.interp_data_path` to a directory with the output
from interpolate_irred.
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse sTDEP parameters.")

    parser.add_argument("--config", type=str, required=True, help="Yaml file containing configuration. See keys in \
                            /algorithms/cv_free_energy HeatCapFreeEnergyParams. Must also contain ucposrcar_path and ssposcar_path.")
    return parser.parse_args()

def main():

    args = parse_arguments()

    with open(args.config) as f:
        cfg_data = yaml.safe_load(f)

    paths = Paths(**cfg_data["Paths"])
    params = HeatCapFreeEnergyParams(**cfg_data["HeatCapFreeEnergyParams"])
    if params.interp_settings is not None:
        params.interp_settings = InterpolateIFCParams(**params.interp_settings)
    else:
        raise ValueError("Expected key interp_settings in HeatCapFreeEnergyParams")

    if params.interp_settings.lds is not None:
        params.interp_settings.lds = LammpsDynamicsSettings(**params.interp_settings.lds)

    if not os.path.isdir(paths.basepath):
        raise RuntimeError("basepath is not a directory")

    os.chdir(paths.basepath)

    setup_logging("sTDEP.log", paths.basepath)
    logging.info(f"CWD: {os.getcwd()}")

    # Move and rename usposcar and ssposcar
    shutil.copyfile(paths.ucposcar_path, os.path.join(paths.basepath, "infile.ucposcar"))
    shutil.copyfile(paths.ssposcar_path, os.path.join(paths.basepath, "infile.ssposcar"))

    run_interpolate = paths.interp_data_path is None
    temp_stencil = initialize_free_energy(params, run_interpolate)

    if run_interpolate:
        _, ss_out_dir = run_interpolate_irred(params.interp_settings, paths)
    else:
        # TODO Check if prior interpolation has the right inputs? 
        # not bothered right now
        logging.info("Using previous interpolation data at {paths.interp_data_path}")
        logging.warning("ASSUMING INTERP WAS DONE WITH PROPER SETTINGS: rc3, rc4 > 0, make_ss_ifces = True, interpolate_U0 = True")

        ss_out_dir = 


    run_cv_free_energy(params, paths, temp_stencil, ss_out_dir)

if __name__ == "__main__":
    main()