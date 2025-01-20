# Ethan's TDEP Tools

The scripts in this repo implement some common algorithms and patterns of using TDEP.

To install the package and its dependencies you can use a virtual environment:
```python
python -m venv venv
source venv/bin/activate
pip install -e .
```

or a conda environment (haven't test this but it should work):
```python
conda create --name <env-name> --file requirements.txt
pip install -e .
```

### Scripts:
- `ifc_from_MD` : Runs a LAMMPS simulation and extracts force constants.
- `sTDEP` : Generates atomic configurations using `canonical_configurations` and self-consistently calculates force constants for a given number iterations. The inital force constants are generated from a maximum frequency. Diagnostics such as the dispersion or density of states are automatically plotted for each iteration to check convergence.
- `interpolate_irred` : Calculates/Loads force constants from a list of temperatures and finds the lagrange/linear interpolant through the irreducible force constants at each temperature. Optionally, the interpolated irreducible force constants can be rebuilt to be comensurate with the original super cell.
- `cv_free_energy` : Calcualtes heat capacity as the second derivative of the `anharmonic_free_energy` with respect to temperature. Force constants are interpolated to each temperature on the finite difference stencil.


Each script is available on the command line when the environment `pip install -e .` was run in is active. To run the scripts simply type:
```
<script-name> --config <path-to-config-yml>
````
There are several example config.yml files in the `./configs` folder and all configurations with documentation (TODO) can be found in `./src/algorithms/configs.py`. Each `config.yml` should have a `Paths` section and another section for the parameters of the command you want to run (e.g. `IFC_MD_Params`).