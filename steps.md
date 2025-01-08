

# Subroutines:


### [Self-Consistent Sampling](#self-consistent-sampling)
    - Inputs: 
        - Super-Cell Size
        - Force Constant Cutoff
        - IFCs OR Max Frequency
        - Tempearture
        - Quantum or Classical
        - Force calculator
    - Outputs : 
        - Finite tempearture configurations
        - Fininte tempearture IFCs
    - Important Checks: 
        - DOS convergence w.r.t. SCF iterations (this only really checks 2nd order conv)
        - IFC convergence w.r.t number of samples (use --stride option)
    - Plots:
        - DOS vs SCF Iteration (Typically < 10 iterations needed)
        - Phonon free energy (just harmonic part???)
    1. Pass initial IFCs or maximum frequency to `canonical_configurations`
    2. Calculate DOS
    3. Repeat until DOS is converged. 

### [Converge Supercell Size and Force Constant Cutoff](#converge-supercell-size-and-force-constant-cutoff)
    - Inputs :
        - Initial supercell size (overestimate ideally)
        - List of cutoffs to test
        - Number of atomic configurations to generate
        - Tempearture
        - Quantum or Classical
        - Force calculator
    - Outputs :
        - Converged cutoff
        - Minimum supercell size required
    - Plots:
        - R^2 vs r-cut

    1. Start with smallest r-cut
    2. Calculate force constants via [sTDEP](#self-consistent-sampling)
    3. Parse R^2 value from cross-validation section log file for second + third + fourth order (ignore polar part)
    4. Repeat 2-3 for all r-cut
    5. Plot of R^2 vs r-cut, choose r-cut that converges R^2 (i.e. captures all interactions)
    6. Choose smallest super-cell size that contains r-cut (r-cut < L/2)

### [Converge K Mesh](#converge-k-mesh)
    - Inputs :
        - 

### [Find Equilibrium Lattice Parameter](#converge-k-mesh)
    DO WE CARE ABOUT DOING THIS? (e.g. accounting for thermal expansion)


### [Determine Finite Difference Scheme](#determine-finite-difference-scheme)
    - Inputs:
        - dT, stencil order
        - IFCs (second, third, fourth)
        - Central temeparture
        - Quantum or Classical
    - Outputs:
        - Plots of bulk heat capacity vs dT for each stencil order

    1. Run unmodified anharmonic free energy to fourth order for all temperatures on all stencils
    2. Calculate heat capacity from $$-T * \frac{\partial ^2F}{\partial T^2}$$

### [Mode Heat Capacity](#mode-heat-capacity)
    - Inputs :
        - dT, stencil order
        - k-mesh resolution
        - IFCs (second, third, fourth)
        - Tempearture
        - Quantum or Classical
    - Outputs :
        - Heat capacity of each phonon mode
    
    1. Run modified anharmonic free energy to fourth order for all temperatures in stencil
    2. Calculate heat capacity from free energy for each mode
        - NEED TO FIGURE OUT HOW TO DISTRIBUTE THE CUMULANT STUFF

# Overall Process

1. Determine super-cell size and force cosntant cutoff at highest temperature (run [this]((#converge-super-cell-size-and-force-constant-cutoff))). We will use the same for all lower temperatures.
    - We can use the sTDEP samples needed for this to determine number of self-consistent iterations to converge. This number can be used in all remaining steps (expect less than 10)
    - We must specify the number of atomic configurations but have not checked convergence of this yet. Choose a large-ish number and verify that it was ok in the next step.
2. Check convergence of force constants with respect to number of samples. Again only test the highest tempearture.
    - Use # SCF iterations as determined by step 1
    - Verify the number of configurations used in step 1 was enough for convergence
    - To check this convergence we should ideally calculate our parameter of interest with the IFCs form each iteration
        - Calculate free energy from harmonic part to verify convergence of the second-order IFCs
        - I propose we check the converngence of U0 = <U - UTEP> which `extract_forceconstants` should dump so we have some idea as to whats going on with the higher order IFCs. U0 is still used in the anharmoinc free energy calculation and is important to converge. Unlike the IFCs which have 3*N_atoms samples, U0 only gets one sample per configuration.
            - This avoids explicitly calculating anharominc free energy to check convergence. Which in principle we should do but at least this gives us an idea. 
3. Find required k-mesh resolution to converge anharmonic free energy calculation (all orders).
4. Determine the effect of finite difference stencil and size of dT.
5. Calculate mode heat capacities using the parameters determined by 1-4




using Symbolics, SymbolicUtils, DynamicDiff, DynamicExpressions
operators = OperatorEnum(; binary_operators=(+, *, /, -), unary_operators=(sin, cos));
variable_names = ["x1", "x2", "x3"];
vars = [Expression(Node{Float64}(feature=i); operators, variable_names) for i in 1:3];

f = vars[1] * sin(vars[2] - 0.5)
deriv_x1 = D(f, 1)
deriv_x1_sym =  convert(SymbolicUtils.Symbolic, deriv_x1, operators; variable_names = variable_names)

myf = build_function(deriv_x1_sym, variable_names...)