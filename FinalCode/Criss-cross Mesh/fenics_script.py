# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 15:37:56 2021

@author: giuli
"""

from dolfin import *
from dolfin import UnitSquareMesh, XDMFFile
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
    
    # Boundary conditions
    [bc.zero(B) for bc in bcs]
    
    # Solver setting
    solver = SLEPcEigenSolver(A, B)
    solver.parameters["solver"] = "krylov-schur"
    solver.parameters["problem_type"] = "gen_hermitian"
    solver.parameters["spectrum"] = "target magnitude"
    solver.parameters["spectral_transform"] = "shift-and-invert"
    
    # Shift Setting
    solver.parameters["spectral_shift"] = 22.0 # 18.5 for P1
    
    # Number of desired eigenvalues to compute
    neigs = 40 
    
    # Solver
    solver.solve(neigs)
    
    # Return the computed eigenvalues in a sorted array
    computed_eigenvalues = []
    for i in range(min(neigs, solver.get_number_converged())):
        r, _ = solver.get_eigenvalue(i) # ignore the imaginary part
        computed_eigenvalues.append(r)
    return np.sort(np.array(computed_eigenvalues))

def print_eigenvalues(mesh):
    
    # Nedelec Elements
    nedelec_V   = FunctionSpace(mesh, "N1curl", 1)
    nedelec_bcs = [DirichletBC(nedelec_V, Constant((0.0, 0.0)), DomainBoundary())]
    nedelec_eig = eigenvalues(nedelec_V, nedelec_bcs)
    
    # P1 Elements
    lagrange_V   = VectorFunctionSpace(mesh, "CG", 1)
    lagrange_bcs = [DirichletBC(lagrange_V.sub(1), 0, "near(x[0], 0) || near(x[0], pi)"),
                    DirichletBC(lagrange_V.sub(0), 0, "near(x[1], 0) || near(x[1], pi)")]
    lagrange_eig = eigenvalues(lagrange_V, lagrange_bcs)
    
    # Exact solution
    true_eig = np.sort(np.array([float(m**2 + n**2) for m in range(7) for n in range(7)]))[1:41]
    
    # Error computation
    errN = np.max(np.abs(nedelec_eig-true_eig))
    errP1 = np.max(np.abs(lagrange_eig-true_eig))
    
    np.set_printoptions(formatter={'float': '{:5.6f}'.format})
    print("Nedelec:  {}".format(nedelec_eig))
    print("P1:       {}".format(lagrange_eig))
    print("Exact:    {}".format(true_eig))
    print("ErrorN:   {}".format(errN))
    print("ErrorP1:  {}".format(errP1))

# Refinement setting    
i = 7
N = 2**i

# Mesh generation
mesh = RectangleMesh(Point(0, 0), Point(pi, pi), N, N, "crossed")

# Mesh storage and export
xdmf = XDMFFile("mesh.xdmf")
xdmf.write(mesh)
xdmf.close()
import meshio
meshio_mesh = meshio.read("mesh.xdmf")
meshio.write("meshio_mesh.msh", meshio_mesh)

print("\ncrossed mesh")
print_eigenvalues(mesh)
    