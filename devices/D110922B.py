from dolfin import *
import numpy as np
import scipy.integrate as integ
import scipy.constants as const
import matplotlib.pyplot as plt
import time, sys, os
my_path = os.path.dirname(os.path.realpath(__file__))
my_path_src = os.path.join(my_path, 'src')

# Import simulation parameter, functions and material constants from following path
sys.path.append(my_path_src)
from sp_functions import *
from material_const import *
from simulation_parameters import *


#### Input parameters from command line ####
fraction_in_dx_centers = float(sys.argv[1])
#fraction_of_free_charges_on_surface = float(sys.argv[2])

## calculate fraction of free charges going to surface with charge neutrality
# desired values in well -2.75376e15
fraction_of_free_charges_on_surface = (-2.75376+(1-fraction_in_dx_centers)*55.53)/((1-fraction_in_dx_centers)*55.53)



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

class Doping_N(Expression):
	def eval(self, value, x):
		value[0] = exp(-(x[0]-d_pos1)*(x[0]-d_pos1)/2.0/(d_width1*d_width1))*dop_1/(d_width1*np.sqrt(2*np.pi)) + \
					exp(-(x[0]-d_pos2)*(x[0]-d_pos2)/2.0/(d_width2*d_width2))*dop_2/(d_width2*np.sqrt(2*np.pi))

class Surface_N(Expression):
	def eval(self, value, x):
		value[0] = -exp(-x[0]/s_width)*surf/s_width



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

ss = dmax/nel*1e-9 				# Step size for integration
gs = dmax/nel					# Grid spacing


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
s_width = .5


print 'grid spacing =', gs
print 'exchange correlation term =', exchange_correlation_term
print 'charge neutral =', charge_neutral
print 'fraction_in_dx_centers = ', fraction_in_dx_centers
print 'fraction_of_free_charges_on_surface = ', fraction_of_free_charges_on_surface
print 'magnitude surface charge = ', surf

