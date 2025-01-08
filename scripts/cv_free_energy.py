

## Steps for single tempearture, T:
# 1. Generate force constants self-consistently with sTDEP
    # Check convergence w.r.t. number of samples, supercell size, cutoff, numer of scf iterations
# 2. Run anharmonic_free_energy at T
    # Check convergence w.r.t k-mesh
# 3. Run anharmonic_free_energy at finite diff points and k-mesh from step 2
 