#### Considered Stucture ####
STRUCTURE = 'D110620C'
#STRUCTURE = 'D110620C_m2'


#### Global simulation parameters ####
nel = 10000							# Number of elements
nomaxit = 1500 					        	# maximum number of iteration
nocs = 3							# number of considered states eigenstates for electron density calculation
t = 10e-3 							# temperature in K
alpha = 0.1	 						# initial mixing factor for poisson solution
interactive_alpha = False 			                # interactively adjusting alpha
beta = 0.1   						        # initial mixing factor for electron densitiy solution
interactive_beta = False 	                		# interactively adjusting beta
target_error_p = 1e-6 				                # traget error for poisson equation in V
target_error_d = 1e18 				                # traget error for electron density in m^-3
exchange_correlation_term = True 	                        # use exchange correlation term
charge_neutral = True				                # if true Fermi-level is calculated considering charge neutraliy
e_fix = 0.							# fix Fermi-level (requires charge_neutral = False)
surface_charge_on = True					# include surface charges in poisson equation


#### Linear Solvers ####
'''
choose linear solver algorithm                            
    "petsc", 'PETSc builtin LU solver'
    "cg", 'Conjugate gradient method',
    "gmres", 'Generalized minimal residual method',
    "minres", 'Minimal residual method',
    "tfqmr", 'Transpose-free quasi-minimal residual method',
    "richardson", 'Richardson method',
    "bicgstab", 'Biconjugate gradient stabilized method'
'''
linear_solver = "petsc"


#### Eigenvalue Solvers ####
'''
choose eigenvalue solver algorithm
    "arnoldi" (Arnoldi) 
    "krylov-schur" (Krylov-Schur) 
'''
slepc_eigensolver = "krylov-schur"


#### Boundary Condition ####
'''
Types of Boundary condition:
        d_vn: Dirichlet at lhs of structure and von Neumann at rhs
        vn_d: Dirichlet at rhs of structure and von Neumann at lhs
'''
BCT = 'd_vn'					                # boundary condition type
elhs = 0						        # energy of conduction band lhs for BCT = 'd_vn' in eV
erhs = 0						        # energy of conduction band rhs for BCT = 'vn_d' in eV


#### Define region where to solve schroedinger equation ####
dqleft = 200                                                    # left boundary of quantum region in nm
dqright = 500                                                   # right boundary of quantum region in nm


#### Debug ####
DEBUG = False                                                   # debug mode. DEBUG = False -> no plots
DEBUG_level = 1					                # Debug level 1: some plots, Debug level 2: all plots
