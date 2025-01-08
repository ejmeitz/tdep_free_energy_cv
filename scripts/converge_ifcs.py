import os
import numpy as np
from typing import List, Optional

def analyze_ifc_cutoffs(
        rc2s : List[float],
        rc3s : Optional[List[float]] = None,
        rc4s : Optional[List[float]] = None,
        n_sTDEP_iters : int = 10,
        n_sTDEP_configs : int = 250
    ) -> float:
    """
    Returns cross validation parameters for the cutoff combinations passed.

    Parameters:
    r_cuts : List[List[float]] : The r_cut values to try, r_cuts[0] gives the list of
        cutoffs to try for second order and so on. Not all combinations will be tested, 
        instead the first set tested will be rcuts[0][0], rcuts[1][0] ...

    Steps:
    1. Initialize r-cut
    2. Calculate force constants via [sTDEP](#self-consistent-sampling)
    3. Parse R^2 value from cross-validation section log file for
        second
        seonc + third
        second + third + fourth order
    4. Repeat 2-3 for all r-cut
    5. Plot of R^2 vs r-cut, choose r-cut that converges R^2 (i.e. captures all interactions)
    6. Choose smallest super-cell size that contains r-cut (r-cut < L/2)
    """
    
    
    all_r_cuts = []
    lens = []
    for rc in [rc2s, rc3s, rc4s]:
        if rc is not None:
            all_r_cuts.append(rc)
            lens.append(len(rc))
        
    # Make sure same number passed for each
    if len(set(lens)) != 1:
        raise ValueError(f"analyze_ifc_cutoffs expects all r-cut lists to have same length. Got : {lens}")
    
    N = len(rc2s) # same as rc3 and rc4, we just checked

    for i in range(N):
        sTDEP()
        parse_cross_validation_output()
        



def converge_ifc_nsamples():
    pass