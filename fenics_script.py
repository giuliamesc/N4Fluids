# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 15:37:56 2021

@author: giuli
"""

from dolfin import *
import numpy as np
if not has_linear_algebra_backend("PETSc"):
    print("DOLFIN has not been configured with PETSc. Exiting.")
    exit()
if not has_slepc():
    print("DOLFIN has not been configured with SLEPc. Exiting.")
    exit()
    
def eigenvalues(V, bcs):
    # Define the bilinear forms on the right- and left-hand sides
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(curl(u), curl(v))*dx
    b = inner(u, v)*dx
    
    # Assemble into PETSc matrices
    dummy = v[0]*dx
    A = PETScMatrix()
    assemble_system(a, dummy, bcs, A_tensor=A)
    B = PETScMatrix()
    assemble_system(b, dummy, bcs, A_tensor=B)
    
    [bc.zero(B) for bc in bcs]
    solver = SLEPcEigenSolver(A, B)
    solver.parameters["solver"] = "krylov-schur"
    solver.parameters["problem_type"] = "gen_hermitian"
    
    solver.parameters["spectrum"] = "target magnitude"
    solver.parameters["spectral_transform"] = "shift-and-invert"
    solver.parameters["spectral_shift"] = 5.5
    neigs = 12
    solver.solve(neigs)
    
     # Return the computed eigenvalues in a sorted array
    computed_eigenvalues = []
    for i in range(min(neigs, solver.get_number_converged())):
        r, _ = solver.get_eigenvalue(i) # ignore the imaginary part
        computed_eigenvalues.append(r)
    return np.sort(np.array(computed_eigenvalues))

def print_eigenvalues(mesh):
    nedelec_V   = FunctionSpace(mesh, "N1curl", 1)
    nedelec_bcs = [DirichletBC(nedelec_V, Constant((0.0, 0.0)), DomainBoundary())]
    nedelec_eig = eigenvalues(nedelec_V, nedelec_bcs)
    
    lagrange_V   = VectorFunctionSpace(mesh, "Lagrange", 1)
    lagrange_bcs = [DirichletBC(lagrange_V.sub(1), 0, "near(x[0], 0) || near(x[0], pi)"),
                    DirichletBC(lagrange_V.sub(0), 0, "near(x[1], 0) || near(x[1], pi)")]
    lagrange_eig = eigenvalues(lagrange_V, lagrange_bcs)
    
    true_eig = np.sort(np.array([float(m**2 + n**2) for m in range(6) for n in range(6)]))[1:13]
    
    np.set_printoptions(formatter={'float': '{:5.2f}'.format})
    print("Nedelec:  {}".format(nedelec_eig))
    print("Lagrange: {}".format(lagrange_eig))
    print("Exact:    {}".format(true_eig))
    
mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40)
print("\ndiagonal mesh")
print_eigenvalues(mesh)

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40, "crossed")
print("\ncrossed mesh")
print_eigenvalues(mesh)
    