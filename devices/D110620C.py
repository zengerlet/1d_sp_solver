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
fraction_of_free_charges_on_surface = float(sys.argv[2])

## calculate fraction of free charges going to surface with charge neutrality
# desired values in well -2.e15
# fraction_of_free_charges_on_surface = (-2.+(1-fraction_in_dx_centers)*12.69)/((1-fraction_in_dx_centers)*12.69)



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


width1 = 10.0       # GaAs
width2 = 310.0      # Al244Ga756As
width3 = 40.0       # GaAs
width4 = 22.0       # Al244Ga756As
width5 = 225.7      # Al2379Ga7621As
width6 = 50.0       # Al1952Ga8048As
width7 = 1000.0     # Al1708Ga8292
wdith8 = 500.0      # GaAs
class Effective_Mass(Expression):
    def eval(self, value, x):
        if between(x[0], layer1):
            value[0] = GaAsff['meff']
        if between(x[0], layer2):
            value[0] = Al244Ga756Asff['meff']
        if between(x[0], layer3):
            value[0] = GaAsff['meff']
        if between(x[0], layer4):
            value[0] = Al244Ga756Asff['meff']
        if between(x[0], layer5):
            value[0] = Al2379Ga7621Asff['meff']
        if between(x[0], layer6):
            value[0] = Al1952Ga8048Asff['meff']
        if between(x[0], layer7):
            value[0] = Al1708Ga8292Asff['meff']
        if between(x[0], layer8):
            value[0] = GaAsff['meff']
        
class Epsilon(Expression):
    def eval(self, value, x):
        if between(x[0], layer1):
            value[0] = GaAsff['eps']
        if between(x[0], layer2):
            value[0] = Al244Ga756Asff['eps']
        if between(x[0], layer3):
            value[0] = GaAsff['eps']
        if between(x[0], layer4):
            value[0] = Al244Ga756Asff['eps']
        if between(x[0], layer5):
            value[0] = Al2379Ga7621Asff['eps']
        if between(x[0], layer6):
            value[0] = Al1952Ga8048Asff['eps']
        if between(x[0], layer7):
            value[0] = Al1708Ga8292Asff['eps']
        if between(x[0], layer8):
            value[0] = GaAsff['eps']


class Band_Energy(Expression):
    def eval(self, value, x):
        if between(x[0], layer1):
            value[0] = GaAsff['cb_e']
        if between(x[0], layer2):
            value[0] = Al244Ga756Asff['cb_e']
        if between(x[0], layer3):
            value[0] = GaAsff['cb_e']
        if between(x[0], layer4):
            value[0] = Al244Ga756Asff['cb_e']
        if between(x[0], layer5):
            value[0] = Al2379Ga7621Asff['cb_e']
        if between(x[0], layer6):
            value[0] = Al1952Ga8048Asff['cb_e']
        if between(x[0], layer7):
            value[0] = Al1708Ga8292Asff['cb_e']
        if between(x[0], layer8):
            value[0] = GaAsff['cb_e']
            

class Doping_N(Expression):
    def eval(self, value, x):
        value[0] = exp(-(x[0]-d_pos1)*(x[0]-d_pos1)/2.0/(d_width1*d_width1))*dop_1/(d_width1*np.sqrt(2*np.pi))

class Surface_N(Expression):
    def eval(self, value, x):
        value[0] = -exp(-x[0]/s_width)*surf/s_width



#### Define Layers ####

width1 = 10.0       # GaAs
width2 = 310.0      # Al244Ga756As
width3 = 40.0       # GaAs
width4 = 22.0       # Al244Ga756As
width5 = 225.7      # Al2379Ga7621As
width6 = 50.0       # Al1952Ga8048As
width7 = 1000.0     # Al1708Ga8292
width8 = 500.0      # GaAs


dmax = width1 + width2 + width3 + width4 + width5 + width6 + width7 + width8
layer1 = (0, width1)
layer2 = (width1, width1+width2)
layer3 = (width1+width2, width1+width2+width3)
layer4 = (width1+width2+width3, width1+width2+width3+width4)
layer5 = (width1+width2+width3+width4, width1+width2+width3+width4+width5)
layer6 = (width1+width2+width3+width4+width5, width1+width2+width3+width4+width5+width6)
layer7 = (width1+width2+width3+width4+width5+width6, width1+width2+width3+width4+width5+width6+width7) 
layer8 = (width1+width2+width3+width4+width5+width6+width7, width1+width2+width3+width4+width5+width6+width7+width8)



ss = dmax/nel*1e-9               # Step size for integration
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

