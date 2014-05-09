#### Considered Stucture ####
STRUCTURE = 'D110922B'


#### Global simulation parameters ####

nel = 40000							# Number of elements
nomaxit = 1500 						# maximum number of iteration
nocs = 3							# number of considered states
t = 1.3 							# temperature
alpha = 0.5	 						# initial mixing factor for poisson solution
interactive_alpha = True 			# interactively adjusting alpha
beta = 1.   						# initial mixing factor for eDensity solution
interactive_beta = False 			# interactively adjusting beta
target_error_p = 1e-6 				# traget error for poisson equation
target_error_d = 1e18 				# traget error for eDensity
exchange_correlation_term = True 	# use exchange correlation term
charge_neutral = True				# if true Fermi-level is calculated considering charge neutraliy
e_fix = 0.							# fix Fermi-level (requires charge_neutral = False)


#### Boundary Condition ####
'''
Types of Boundary condition:
	d_vn: Dirichlet at lhs of structure and von Neumann at rhs
	vn_d: Dirichlet at rhs of structure and von Neumann at lhs
'''
BCT = 'd_vn'					# boundary condition type
elhs = 0						# energy of conduction band lhs for BCT = 'd_vn'
erhs = 0						# energy of conduction band rhs for BCT = 'vn_d'


#### Define region where to solve schroedinger equation ####
dqleft = 200
dqright = 280

#### Debug ####

DEBUG = True					# debug mode. DEBUG = False -> no plots
DEBUG_level = 1					# Debug level 1: some plots, Debug level 2: all plots
