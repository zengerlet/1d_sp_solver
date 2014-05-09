#### Global simulation parameters ####

nel = 10000						# Number of elements
nomaxit = 1500 					# maximum number of iteration
nocs = 3						# number of considered states
t = 1.3 						# temperature
v_schottky = 0.53				# schottky barrier
alpha = 0.5	 					# initial mixing factor for poisson solution
interactive_alpha = True 		# interactively adjusting alpha
beta = 1.   					# initial mixing factor for eDensity solution
interactive_beta = False 		# interactively adjusting beta
target_error_p = 1e-6 			# traget error for poisson equation
target_error_d = 1e18 			# traget error for eDensity
# Define region where to solve schroedinger equation
dqleft = 200
dqright = 280
exchange_correlation_term = False 	# use exchange correlation term
charge_neutral = True				# if true Fermi-level is calculated considering charge neutraliy
e_fix = 0.							# fix Fermi-level (requires charge_neutral = False)

DEBUG = False					# debug mode. DEBUG = False -> no plots
DEBUG_level = 1
