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
# desired values in well -2.e15
fraction_of_free_charges_on_surface = (-2.+(1-fraction_in_dx_centers)*12.69)/((1-fraction_in_dx_centers)*12.69)



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
            value[0] = GaAsnn['meff']
        if between(x[0], layer4):
            value[0] = Al24Ga76Asnn['meff']
        
class Epsilon(Expression):
    def eval(self, value, x):
        if between(x[0], layer1):
            value[0] = GaAsnn['eps']
        if between(x[0], layer2):
            value[0] = Al24Ga76Asnn['eps']
        if between(x[0], layer3):
            value[0] = GaAsnn['eps']
        if between(x[0], layer4):
            value[0] = Al24Ga76Asnn['eps']


class Band_Energy(Expression):
    def eval(self, value, x):
        if between(x[0], layer1):
            value[0] = GaAsnn['cb_e']
        if between(x[0], layer2):
            value[0] = Al24Ga76Asnn['cb_e']
        if between(x[0], layer3):
            value[0] = GaAsnn['cb_e']
        if between(x[0], layer4):
            value[0] = Al24Ga76Asnn['cb_e']
            

class Doping_N(Expression):
    def eval(self, value, x):
        value[0] = exp(-(x[0]-d_pos1)*(x[0]-d_pos1)/2.0/(d_width1*d_width1))*dop_1/(d_width1*np.sqrt(2*np.pi))

class Surface_N(Expression):
    def eval(self, value, x):
        value[0] = -exp(-x[0]/s_width)*surf/s_width



#### Define Layers ####

dmax = 1625.7

layer1 = (0.0,            # GaAs
          10.)    
layer2 = (10.,         # AlGaAs
          320.)    
layer3 = (320.,             # GaAs
          360.)            
layer4 = (360.,             #AlGaAs
          dmax)    



ss = dmax/nel*1e-9                 # Step size for integration
gs = dmax/nel                    # Grid spacing


#### Doping ####
# Define doping layers
d_pos1 = 250.
d_width1 = .15

#fraction_in_dx_centers = 0.892
#fraction_in_dx_centers = 0
dop_1 = (1.0-fraction_in_dx_centers)*12.69e24                 


### Surface charge ###
#fraction_of_free_charges_on_surface = 0.8
surf = fraction_of_free_charges_on_surface*(dop_1)
s_width = .5


print 'grid spacing =', gs
print 'exchange correlation term =', exchange_correlation_term
print 'charge neutral =', charge_neutral
print 'fraction_in_dx_centers = ', fraction_in_dx_centers
print 'fraction_of_free_charges_on_surface = ', fraction_of_free_charges_on_surface

