from dolfin import *
import numpy as np
import scipy.integrate as integ
import scipy.constants as const
import matplotlib.pyplot as plt
import time, sys, os
my_path = os.path.dirname(os.path.realpath(__file__))
my_path_src = os.path.join(my_path, 'src')
my_path_devices = os.path.join(my_path, 'devices')

# Parameter for Fenics
parameters['reorder_dofs_serial'] = False

# Import simulation parameter, functions and material constants from following path
sys.path.append(my_path_src)
from sp_functions import *
from material_const import *
from simulation_parameters import *

# Import structure
sys.path.append(my_path_devices)
if STRUCTURE == 'D110922B':
    from D110922B import *
if STRUCTURE == 'STRUCT1':
    from STRUCT1 import *
if STRUCTURE == 'D110620C':
    from D110620C import *
if STRUCTURE == 'D110620C_m2':
    from D110620C_m2 import *


# get time program runs
start_time = time.time()

#### Built framework for FEM solver ####

# Initialize sub-domain instances
boundary_left = Boundary_Left()
boundary_right = Boundary_Right()
q_Domain = Q_Domain()

# Define mesh
mesh = IntervalMesh(nel,0.0,dmax)

# array for mesh
x = mesh.coordinates()
np.savetxt('x_array_' + str(nel) + '.out',x)


# Initialize mesh function for boundaries and domains
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
boundary_left.mark(boundaries, 1)
boundary_right.mark(boundaries, 2)
domains = CellFunction("size_t", mesh)
domains.set_all(0)
q_Domain.mark(domains, 1)

# Initialize submesh for quantum region
q_mesh = SubMesh(mesh, domains, 1)
q_domain = CellFunction("size_t", q_mesh)
x_q = q_mesh.coordinates()
np.savetxt('x_q_array_' + str(nel) + '.out',x_q)

nelq = np.size(x_q) - 1        # number of elements in quantum region

# Define function space and basis function for schroedinger equation
Vq = FunctionSpace(q_mesh, "Lagrange", 1)
uq = TrialFunction(Vq)
vq = TestFunction(Vq)
pot_tot_q = Function(Vq)
eDens_q = Function(Vq)

# Define function space and basis functions for poisson equation
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)

# Define constants and functions for variational for of Schroedinger equation
h1 = Constant(const.hbar**2/2/const.m_e/const.e*1e18)  # 1e18 is to convert m into nm
m_eff = Effective_Mass()
m_eff_array_q = project(m_eff, Vq).vector().array()
m_eff_array = project(m_eff, V).vector().array()
bandE = Band_Energy()
bandE_array_q = project(bandE, Vq).vector().array()
bandE_array = project(bandE, V).vector().array()
np.savetxt('bandE_array_' + str(nel) + '.out',bandE_array)

if DEBUG and (DEBUG_level==1 or DEBUG_level==2):
    plt.figure()
    plt.plot(x, bandE_array)
    plt.title('Conduction band')


# Define constants and functions for variational form of Poisson equation
c1 = Constant('1e18') # convert m into nm
c2 = Constant(const.e/const.epsilon_0)
epsilon = Epsilon()
epsilon_array_q = project(epsilon, Vq).vector().array()
epsilon_array = project(epsilon, V).vector().array()
np.savetxt('epsilon_array_' + str(nel) + '.out',epsilon_array)
doping_n = Doping_N()
doping_n_array = project(doping_n, V).vector().array()
print 'doping net charge', netCharge(doping_n_array,ss)
np.savetxt('doping_n_array_' + str(nel) + '.out',doping_n_array*1e-24)
surface_charge = Surface_N()
surface_charge_array = project(surface_charge, V).vector().array()
print 'surface net charge', netCharge(surface_charge_array,ss)

##################################


if DEBUG and (DEBUG_level==1 or DEBUG_level==2):
    plt.figure()
    plt.plot(x, doping_n_array)
    plt.title('Doping')
    plt.figure()
    plt.plot(x, surface_charge_array)
    plt.title('Surface charge')
    plt.show()


# Define Dirichlet boundary conditions for Poisson equation at left and right boundaries
if BCT == 'd_vn':
    bcs = [DirichletBC(V, 0.0, boundaries, 1)]

if BCT == 'vn_d':
    bcs = [DirichletBC(V, 0.0, boundaries, 2)]

# slope on rhs
g_R = Constant('0.0')
   
# Define measures
ds = Measure("ds")[boundaries]
dx = Measure("dx")[domains]
dx_q = Measure("dx")[q_domain]

# functions for solutions
eDens1 = Function(Vq)
eDens2 = Function(Vq)
eDensp = Function(V)
u_p1 = Function(V)
u_p2 = Function(V)
pot_tot = Function(Vq)

# initialize array for solution of poisson equation
u_p1_array = np.zeros(nel+1)
# initialize array for solution of electron density
eDens1_array = np.zeros(nelq+1)
# initialize array for solution of schroedinger equation
Psi = np.zeros([nocs, nelq+1])


#### Schroedinger Poisson ####
noit = 0
error_p1 = 0
error_p2 = 0
ierror_p = 0
error_d1 = 0
error_d2 = 0
ierror_d = 0
while noit < nomaxit:
    
    ###################### Schroedinger ######################
    
    # update pot_tot 
    pot_tot_array_p = bandE_array - u_p1_array
    
    if BCT == 'd_vn':
        # for charge neutral. set lhs of conduction band to 0ev
        pot_tot_array_p = pot_tot_array_p - pot_tot_array_p[0] + elhs
    
    if BCT == 'vn_d':
        # for charge neutral. set lhs of conduction band to 0ev
        pot_tot_array_p = pot_tot_array_p - pot_tot_array_p[-1] + erhs
    
    # go from larger grid from poisson equation to smaller grid for schroedinger equation and include exchange correlation term
    
    if exchange_correlation_term:
        pot_tot_array = pot_tot_array_p[q_mesh.data().array('parent_vertex_indices', 0)] + V_ex(eDens1_array, epsilon_array_q, m_eff_array_q)/const.e
    else: pot_tot_array = pot_tot_array_p[q_mesh.data().array('parent_vertex_indices', 0)]
    
    if DEBUG and DEBUG_level==2:
        # plot u_p1
        plt.figure()
        plt.plot(x, u_p1_array)
        plt.title('from poisson')
        # plot new total potential
        plt.figure()
        plt.plot(x_q, pot_tot_array)
        plt.title('new total potential for schroedinger')
        plt.show()
    
    pot_tot.vector()[:] = np.array(pot_tot_array)

    # Define variational form
    h = (h1*inner(1.0/m_eff*grad(uq), grad(vq))*dx_q(0) + (pot_tot)*uq*vq*dx_q(0))
    m = (uq*vq*dx_q(0))

    # Assemble stiffness form
    H = PETScMatrix()
    M = PETScMatrix()
    assemble(h, tensor=H)
    assemble(m, tensor=M)

    # Create eigensolver
    eigensolver = SLEPcEigenSolver(H,M)
    eigensolver.parameters["spectrum"]="smallest real"
    #eigensolver.parameters["solver"]="lapack"

    # Compute all eigenvalues of A x = \lambda x
    print "Computing eigenvalues. This can take a minute."
    eigensolver.solve(nocs)

    # Print and display all bound states
    r = 0
    dex = 0
    E = []

    while dex < nocs:
        # Extract next smallest eigenpair
        r, c, rx, cx = eigensolver.get_eigenpair(dex)
        Psi1 = normalize(rx, ss)
        Psi[dex,:] = Psi1
        E.append(r)
        dex += 1
    

    print 'eigenvalues =', E
    
    ### Density calculation
    if charge_neutral:
        ef = findFermi(netCharge(doping_n_array, ss)+netCharge(surface_charge_array, ss), Psi, E, m_eff_array_q, t, ss, nelq, nocs)
    if not(charge_neutral):
        ef = e_fix
    eDens2_array = eDensity(Psi, m_eff_array_q, E, ef, t, nelq, nocs, True)
    error_d2 = max(abs(eDens2_array - eDens1_array))
    print 'eDensity error in step =', noit, error_d2
    print 'eDensity target error =', target_error_d
    print 'beta =', beta
    
    # update error
    error_d1 = error_d2

    if DEBUG and DEBUG_level==2:
        for k in xrange(nocs):
            print "E[" + str(k) + "]=" + str(E[k])
        plt.figure()
        plt.plot(x_q, eDens2_array)
        plt.show()

    
    # mix old and new solution with mixing factor alpha (beta = fraction of new solution)
    eDens1_array = beta*eDens2_array + (1.0 - beta)*eDens1_array

    # use dens_array as coefficients function space V
    xind = q_mesh.data().array('parent_vertex_indices', 0)
    eDensp.vector()[xind] = np.array(eDens1_array)
    
    if DEBUG and DEBUG_level==2:
        test_array = project(eDensp, V).vector().array()
        plt.figure()
        plt.plot(x_q,eDens1_array)
        plt.plot(x,test_array)
        plt.show()

    
    ###################### Poisson ######################
        
    # Define variational form
    F = ((doping_n - eDensp + surface_charge)*c2*v*dx(0) + (doping_n - eDensp + surface_charge)*c2*v*dx(1))
    a = (c1*inner(epsilon*grad(u), grad(v))*dx(0) + c1*inner(epsilon*grad(u), grad(v))*dx(1))
    
    # Solve problem
    '''
    solve(a == F, u_p2, bcs)
    '''
    solve(a == F, u_p2, bcs=bcs,
        solver_parameters={"linear_solver": "cg"})    

    # convert solution to numpy array
    u_p2_array = u_p2.vector().array()
    error_p2 = max(abs(u_p2_array - u_p1_array))
    print 'poisson error in step =', noit, error_p2
    print 'poisson target error =', target_error_p
    print 'alpha =', alpha
    
    # interactive alpha
    if interactive_alpha:
        # refine alpha if convergence is not monotonous
        if (error_p2 - error_p1) > 0 and noit != 0:
            alpha = alpha/2
            ierror_p = 0
        
        # make alpha larger if not refined for 5 iterations
        if ierror_p == 5 and alpha*1.2<1:
            alpha = alpha*1.2
            ierror_p = 0

        if (noit == nomaxit-1):
            print 'not converged, final error =',  error_p2
        
        ierror_p += 1
    
    # interactive alpha
    if interactive_beta:
        # refine beta if convergence is not monotonous
        if (error_d2 - error_d1) > 0 and noit != 0:
            beta = beta/2
            ierror_d = 0
        
        # make alpha larger if not refined for 5 iterations
        if ierror_d == 5 and beta*1.2<1:
            beta = beta*1.2
            ierror_d = 0

        if (noit == nomaxit-1):
            print 'not converged, final error =',  error_d2
        
        ierror_d += 1
    
    
    # break if error gets below target error
    if error_p2 < target_error_p and error_p2 != 0 and error_d2 < target_error_d and error_d2 != 0:
        print 'Converged!'
        # stop time
        time_ex = time.time() - start_time
        plot_output(x, x_q, pot_tot_array_p, doping_n_array, eDens1_array, nel, ef, time_ex, 
                    noit, target_error_p, error_p2, nocs, E, Psi, ss, gs, nomaxit, exchange_correlation_term, 
                    DEBUG, DEBUG_level, fraction_in_dx_centers, fraction_of_free_charges_on_surface, surface_charge_array, t)
        break

    # update error
    error_p1 = error_p2
    
    # mix old and new solution with mixing factor alpha (alpha = fraction of new solution)
    u_p1_array = alpha*u_p2_array + (1.0 - alpha)*u_p1_array
    noit += 1


