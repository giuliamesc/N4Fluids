# Eigenvalues of the Maxwell Cavity Problem :zap:
This folder contains all the code and results related to the project 2 for the course MATH-468 Numerics for Fluids, Structures and Electromagnetics at EPFL.

## Triangular Mesh
This folder contains the FreeFem++ codes to solve the problem on a triangular mesh using P1 and Nedelec elements, as well as files containing the computed eigenvalues for both methods.

## Unstructured Mesh
This folder contains the FreeFem++ codes to solve the problem on an unstructured mesh using P1 and Nedelec elements, as well as files containing the computed eigenvalues for both methods.

## Criss-cross Mesh
This folder contains the FeNiCs code to generate the criss-cross mesh. The problem is solved as well for P1 and Nedelec elements. The code prints the computed eigenvalues on the terminal, which have been manually copied in the matlab file `Plots/plot_P1.m ` to plot the results. This code also saves the mesh in the file `meshio_mesh.msh`, which can be visualized with GMSH, and exported in the *.mesh* format. This format can be read by FreeFem++, but first it has to be preprocessed to remove the last column, due to the dimension of the mesh.

## Convergence
This folder contains the code used to run convergence tests on the Unstructured mesh with Nedelec elements and find a suitable refinement level for the next computations. The file `unstructured_Nedelec_convergence.edp` solves the problem with different levels of refinements (variable ***I***) and saves the results in the folder `ConvergenceTests/`, under the name `lambda_unstructured_Nedelec_REFINEMENT.dat`.

## Plots
This folder contains the Matlab files used to plot the results:
* `convtest_unstr_nedelec.m ` plots the results of the convergence test desc.
* `plot_P1.m ` and ` plot_Ned.m ` are used to plot the computed eigenvalues for the best methods using P1 and Nedelec elements against the real ones.

## Authors

* [Giulia Mescolini](https://github.com/giuliamesc) :ghost:
* [Thomas Rimbot](https://github.com/Thomas-debug-creator) 🕷️
