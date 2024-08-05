# exact_diag
Exact diagonalization code in Julia, used to calculate XAS of d-block atoms

This code describes basis states using integers, where the 1s and 0s represent occupied and unoccupied orbitals in (ml, ms) basis.

"params.jl" defines some useful calculation parameters

"tables.jl" defines useful constants, such as Gaunt coefficients and crystal-field matrices

"tools.jl" defines numerous functions that are used for calculating XAS and Tunabe-Sugano curves, which are fairly general and can likely be adapted for further exact diagonalization problems

"XAS.ipynb" is a Jupyter Notebook that shows how the code can be used to calcuate Tunabe-Sugano curves and XAS spectra. 
