from dolfin import *
import numpy as np
import scipy.integrate as integ
import scipy.constants as const
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('/cluster/home03/phys/thomazen/D110922B/src')
from sp_functions import *
from material_const import *
parameters['reorder_dofs_serial'] = False



def plot_output(x, x_q, pot_tot_array_p, doping_n_array, eDens_array, nel, ef, time_ex, noit, target_error, error, nocs, E, PSI, ss, gs, nomaxit, exchange_correlation_term):
	'Output all relevant information and write into file'
	if DEBUG and (DEBUG_level==1 or DEBUG_level==2):
		plt.figure()
		plt.plot(x_q, pot_tot_array)
		for k in xrange(nocs):
			print "E[" + str(k) + "]=" + str(E[k])
			plt.plot(x_q, PSI[k]/np.max(abs(PSI[k])) + E[k])
		plt.title('Conduction band, wave functions, fermi level')	
		plt.show()
		
	if exchange_correlation_term:
		myfile = open('output_ex_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.txt', 'w')
	else: myfile = open('output_'  + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.txt', 'w')
	myfile.write('doping net charge = ' + str(netCharge(doping_n_array, ss)) + '\n')
	myfile.write('eDensity net charge = ' + str(-netCharge(eDens_array, ss)) + '\n') 
	myfile.write('surface net charge = ' + str(netCharge(surface_charge_array, ss)) + '\n') 
	myfile.write('total net charge = ' + str(netCharge(doping_n_array, ss)-netCharge(eDens_array, ss)+netCharge(surface_charge_array, ss)) + '\n')
	myfile.write('fraction_in_dx_centers = ' + str(fraction_in_dx_centers) + '\n')
	myfile.write('fraction_of_free_charges_on_surface = ' + str(fraction_of_free_charges_on_surface) + '\n')
	myfile.write('eigenvalues = ' + str(E) + '\n')
	myfile.write('converged at step ' + str(noit) + '\n')
	myfile.write('max number of iterations = ' + str(nomaxit) + '\n')
	myfile.write('target error = ' + str(target_error) + '\n')
	myfile.write('error = ' + str(error) + '\n')
	myfile.write('number of elements = ' + str(nel) + '\n')
	myfile.write('number of considered states of schroedinger equation = ' + str(nocs) + '\n')
	myfile.write('exchange correlation = ' + str(exchange_correlation_term) + '\n')
	myfile.write('calculation time = ' + str(time_ex) + 'sec\n')
	myfile.write('fermi level = ' + str(ef) + '\n')
	myfile.write('grid spacing = ' + str(gs))
	myfile.close()
	if exchange_correlation_term:
		np.savetxt('eDens_array_ex_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', eDens_array*1e-24)
	else: np.savetxt('eDens_array_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', eDens_array*1e-24)
	if exchange_correlation_term:
		np.savetxt('pot_tot_array_p_ex_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', pot_tot_array_p)
	else: np.savetxt('pot_tot_array_p_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', pot_tot_array_p)
	if exchange_correlation_term:
		np.savetxt('Psi_ex_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', np.hstack([x_q,PSI.T]))
	else: np.savetxt('Psi_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(fraction_of_free_charges_on_surface) + '.out', np.hstack([x_q,PSI.T]))
	print 'doping net charge = ', netCharge(doping_n_array, ss)
	print 'eDensity net charge = ', netCharge(eDens_array, ss)
	print 'total net charge = ', netCharge(doping_n_array, ss)-netCharge(eDens_array, ss)+netCharge(surface_charge_array, ss)
	print 'fraction_in_dx_centers = ', fraction_in_dx_centers
	print 'fraction_of_free_charges_on_surface = ', fraction_of_free_charges_on_surface
	print 'converged at step', noit
	print 'target error =', target_error
	print 'error =', error
	print 'number of elements =', nel
	print 'number of considered states of schroedinger equation =', nocs
	print 'exchange correlation =', exchange_correlation_term
	print 'calculation time =', time_ex
	

# Define sub-domains and boundaries
class Boundary_Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class Boundary_Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], dmax)

class Q_Domain(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (dqleft, dqright))

class Effective_Mass(Expression):
    def eval(self, value, x):
		if between(x[0], layer1):
			value[0] = GaAsnn['meff']
		if between(x[0], layer2):
			value[0] = Al24Ga76Asnn['meff']
		if between(x[0], layer3):
			value[0] = AlAsnn['meff']
		if between(x[0], layer4):
			value[0] = GaAsnn['meff']
		if between(x[0], layer5):
			value[0] = AlAsnn['meff']
		if between(x[0], layer6):
			value[0] = GaAsnn['meff']
		if between(x[0], layer7):
			value[0] = AlAsnn['meff']
		if between(x[0], layer8):
			value[0] = GaAsnn['meff']
		if between(x[0], layer9):
			value[0] = AlAsnn['meff']
		if between(x[0], layer10):
			value[0] = Al24Ga76Asnn['meff']
		if between(x[0], layer11):
			value[0] = GaAsnn['meff']
		if between(x[0], layer12):
			value[0] = Al24Ga76Asnn['meff']
		if between(x[0], layer13):
			value[0] = AlAsnn['meff']
		if between(x[0], layer14):
			value[0] = GaAsnn['meff']
		if between(x[0], layer15):
			value[0] = AlAsnn['meff']
		if between(x[0], layer16):
			value[0] = GaAsnn['meff']
		if between(x[0], layer17):
			value[0] = AlAsnn['meff']
		if between(x[0], layer18):
			value[0] = GaAsnn['meff']
		if between(x[0], layer19):
			value[0] = AlAsnn['meff']
		if between(x[0], layer20):
			value[0] = Al24Ga76Asnn['meff']
		
class Epsilon(Expression):
    def eval(self, value, x):
		if between(x[0], layer1):
			value[0] = GaAsnn['eps']
		if between(x[0], layer2):
			value[0] = Al24Ga76Asnn['eps']
		if between(x[0], layer3):
			value[0] = AlAsnn['eps']
		if between(x[0], layer4):
			value[0] = GaAsnn['eps']
		if between(x[0], layer5):
			value[0] = AlAsnn['eps']
		if between(x[0], layer6):
			value[0] = GaAsnn['eps']
		if between(x[0], layer7):
			value[0] = AlAsnn['eps']
		if between(x[0], layer8):
			value[0] = GaAsnn['eps']
		if between(x[0], layer9):
			value[0] = AlAsnn['eps']
		if between(x[0], layer10):
			value[0] = Al24Ga76Asnn['eps']
		if between(x[0], layer11):
			value[0] = GaAsnn['eps']
		if between(x[0], layer12):
			value[0] = Al24Ga76Asnn['eps']
		if between(x[0], layer13):
			value[0] = AlAsnn['eps']
		if between(x[0], layer14):
			value[0] = GaAsnn['eps']
		if between(x[0], layer15):
			value[0] = AlAsnn['eps']
		if between(x[0], layer16):
			value[0] = GaAsnn['eps']
		if between(x[0], layer17):
			value[0] = AlAsnn['eps']
		if between(x[0], layer18):
			value[0] = GaAsnn['eps']
		if between(x[0], layer19):
			value[0] = AlAsnn['eps']
		if between(x[0], layer20):
			value[0] = Al24Ga76Asnn['eps']

class Band_Energy(Expression):
    def eval(self, value, x):
		if between(x[0], layer1):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer2):
			value[0] = Al24Ga76Asnn['cb_e']
		if between(x[0], layer3):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer4):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer5):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer6):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer7):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer8):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer9):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer10):
			value[0] = Al24Ga76Asnn['cb_e']
		if between(x[0], layer11):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer12):
			value[0] = Al24Ga76Asnn['cb_e']
		if between(x[0], layer13):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer14):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer15):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer16):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer17):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer18):
			value[0] = GaAsnn['cb_e']
		if between(x[0], layer19):
			value[0] = AlAsnn['cb_e']
		if between(x[0], layer20):
			value[0] = Al24Ga76Asnn['cb_e']


'''
class Doping_N(Expression):
    def eval(self, value, x):
        value[0] = exp(-(x[0]-d_pos1)*(x[0]-d_pos1)/2.0/(d_width1*d_width1))*dop_1 + \
					exp(-(x[0]-d_pos2)*(x[0]-d_pos2)/2.0/(d_width2*d_width2))*dop_2

'''


class Doping_N(Expression):
    def eval(self, value, x):
        value[0] = exp(-(x[0]-d_pos1)*(x[0]-d_pos1)/2.0/(d_width1*d_width1))*dop_1/(d_width1*np.sqrt(2*np.pi)) + \
					exp(-(x[0]-d_pos2)*(x[0]-d_pos2)/2.0/(d_width2*d_width2))*dop_2/(d_width2*np.sqrt(2*np.pi))


class Surface_N(Expression):
    def eval(self, value, x):
        value[0] = -exp(-x[0]/s_width)*surf/s_width



# get time program runs
start_time = time.time()


#### Define Layers ####

DopingWidth = 0.5
CapWidth = 10.
OuterSpacerWidth = 133.6
AlAsWidth = 1.981
DopingGaAsWidth = 1.415
InnerSpacerWidth = 69.1
QWellWidth = 30.
DownInnerSpacerWidth = 74.1
DownOuterSpacerWidth = 1738.
#DownOuterSpacerWidth = 500.

layer1 = (0.0,			# GaAs
	      CapWidth)	
layer2 = (CapWidth, 		# AlGaAs
		  CapWidth+OuterSpacerWidth)	
layer3 = (CapWidth+OuterSpacerWidth, 			# AlAs
		  CapWidth+OuterSpacerWidth+AlAsWidth)			
layer4 = (CapWidth+OuterSpacerWidth+AlAsWidth, 			#GaAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth)	
layer5 = (CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth,			# AlAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth)	
layer6 = (CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth, 			# GaAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth)	
layer7 = (CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth, 			# AlAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth)	
layer8 = (CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth, 			# GaAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth)	 
layer9 = (CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth, 			# AlAs
		  CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth)	

AlAs4LowerCoord = CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth

layer10 = (AlAs4LowerCoord, 	# AlGaAs
		   AlAs4LowerCoord+InnerSpacerWidth)	
layer11 = (AlAs4LowerCoord+InnerSpacerWidth,	# GaAs
		   AlAs4LowerCoord+InnerSpacerWidth+QWellWidth)				
layer12 = (AlAs4LowerCoord+InnerSpacerWidth+QWellWidth,		# AlGaAs	
		   AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth)
layer13 = (AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth,	# AlAs
		   AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth)
layer14 = (AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth,		# GaAs
		   AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth+DopingGaAsWidth)
layer15 = (AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth+DopingGaAsWidth,		# AlAs
		   AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth)

AlAs6LowerCoord = AlAs4LowerCoord+InnerSpacerWidth+QWellWidth+DownInnerSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth

layer16 = (AlAs6LowerCoord,		# GaAs
		   AlAs6LowerCoord+DopingGaAsWidth)
layer17 = (AlAs6LowerCoord+DopingGaAsWidth,		# AlAs
		   AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth)
layer18 = (AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth,		# GaAs
		   AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth)
layer19 = (AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth,		# AlAs
		   AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth)
layer20 = (AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth,		# AlGaAs
		   AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DownOuterSpacerWidth)



dmax = AlAs6LowerCoord+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DownOuterSpacerWidth



#### Input parameters from command line ####
fraction_in_dx_centers = float(sys.argv[1])
fraction_of_free_charges_on_surface = float(sys.argv[2])

## calculate fraction of free charges going to surface with charge neutrality
# desired values in well -2.75376e15
fraction_of_free_charges_on_surface = (-2.75376+(1-fraction_in_dx_centers)*55.53)/((1-fraction_in_dx_centers)*55.53)


print 'fraction_in_dx_centers = ', fraction_in_dx_centers
print 'fraction_of_free_charges_on_surface = ', fraction_of_free_charges_on_surface

#### General Parameter ####

#dmax = 450.0 					# Total thickess of structure in nm
nel = 30000						# Number of elements
nomaxit = 1500 					# maximum number of iteration
ss = dmax/nel*1e-9 				# Step size for integration
gs = dmax/nel					# Grid spacing
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
#dqleft = 200
#dqright = 280
dqleft = 100
dqright = 400
exchange_correlation_term = False 	# use exchange correlation term
charge_neutral = True				# if true Fermi-level is calculated considering charge neutraliy
e_fix = 0.							# fix Fermi-level (requires charge_neutral = False)

DEBUG = False					# debug mode. DEBUG = False -> no plots
DEBUG_level = 1


print 'grid spacing =', gs
print 'exchange correlation term =', exchange_correlation_term
print 'charge neutral =', charge_neutral


#### Doping ####

# Define doping layers
d_pos1 = CapWidth+OuterSpacerWidth+AlAsWidth+DopingGaAsWidth+AlAsWidth+DopingGaAsWidth-DopingWidth/2-0.6
d_width1 = .15
d_pos2 = AlAs6LowerCoord+DopingGaAsWidth-DopingWidth/2-0.6
d_width2 = .15
#fraction_in_dx_centers = 0.892
#fraction_in_dx_centers = 0
#dop_1 = (1.0-fraction_in_dx_centers)*86.38e24 				# doping in m^-3 
dop_1 = (1.0-fraction_in_dx_centers)*43.19e24 				
#dop_2 = (1.0-fraction_in_dx_centers)*24.68e24
dop_2 = (1.0-fraction_in_dx_centers)*12.34e24


### Surface charge ###
#fraction_of_free_charges_on_surface = 0.8
surf = fraction_of_free_charges_on_surface*(dop_1+dop_2)
print 'magnitude surface charge = ', surf
s_width = .5


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

nelq = np.size(x_q) - 1		# number of elements in quantum region

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
np.savetxt('doping_n_array_' + str(nel) + '.out',doping_n_array*1e-24)
surface_charge = Surface_N()
surface_charge_array = project(surface_charge, V).vector().array()


print 'next nano net doping = ', 0.606150e12*10000
print 'doping net charge = ', netCharge(doping_n_array, ss)
print 'surface net charge = ', netCharge(surface_charge_array, ss)
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
bcs = [DirichletBC(V, 0.0, boundaries, 1)]

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
	
	'''
	# calculate offset to keep band at left hand side to 0.53eV
	pot_tot_array_p = pot_tot_array_p - pot_tot_array_p[0] + v_schottky
	'''
	# for charge neutral. set lhs of conduction band to 0ev
	pot_tot_array_p = pot_tot_array_p - pot_tot_array_p[0]
	
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
					noit, target_error_p, error_p2, nocs, E, Psi, ss, gs, nomaxit, exchange_correlation_term)
		break

	# update error
	error_p1 = error_p2
	
	# mix old and new solution with mixing factor alpha (alpha = fraction of new solution)
	u_p1_array = alpha*u_p2_array + (1.0 - alpha)*u_p1_array
	noit += 1


