import numpy as np
import scipy.integrate as integ
import scipy.constants as const
import matplotlib.pyplot as plt

def wf_weights(WF_array, EV, nocs, ef, temp):
    """ Helper function to calculate the exponential part of the delectron density """
    WF_weights = np.zeros(nocs)
    # sum up small values (exponential factor for higher energies) first to prevent from numerical errors
    for i in xrange(nocs-1,-1,-1):  
        exparg = -(EV[i] - ef)*const.e/const.k/temp
        # approximate for large arguments
        if exparg > 100.:
            WF_weights[i] = exparg
        else: WF_weights[i] = np.log(1 + np.exp(exparg))
    return WF_weights

def sigma_z(x_array, eDens_array):
    """ Calculation of sigma of electron density """
    x_array = x_array.T[0,:]
    z1 = (integ.simps(x_array*eDens_array,x_array)/integ.simps(eDens_array,x_array))**2
    z2 = integ.simps(x_array**2*eDens_array,x_array)/integ.simps(eDens_array,x_array)
    return np.sqrt(z2-z1)

def normalize(f, ss):
    """ Normalize wave function """
    N = integ.simps(f*f, dx=ss)
    c = 1.0/np.sqrt(N)
    f_norm = c*f
    return f_norm

def netCharge(DOP_array, ss):
    """ Integrate charge density array in z-direction """
    DOP_net_charge = integ.simps(DOP_array, dx=ss) 
    return DOP_net_charge

def eDensity(WF_array, m_eff_array_q, EV, ef, temp, nelq, nocs, oW=False):
    """ Calculate eDensity from wave functions considering fermi-dirac distribution and density of states """
    density = np.zeros(nelq + 1)
    Weights = wf_weights(WF_array, EV, nocs, ef, temp)
    for i in xrange(nocs-1,-1,-1):
        density = density + WF_array[i]**2*Weights[i]
    if oW:
        print 'Weights of wave functions = ', Weights
    # multiply by constants and density of states 
    DOS = m_eff_array_q*const.m_e/const.pi/const.hbar**2
    density = density*temp*const.k*DOS
    return density
    
def findFermi(net_Doping, WF_array, EV, m_eff_array_q, temp, ss, nelq, nocs):
    """ Find fermi level using charge neutrality condition """
    intervalok = True
    #start with ef_min = -1.0, ef_max = 1.0
    ef_min = -1.0
    ef_max = 1.0
    # check if given interval is ok
    D_ef_min = eDensity(WF_array, m_eff_array_q, EV, ef_min, temp, nelq, nocs)
    D_ef_max = eDensity(WF_array, m_eff_array_q, EV, ef_max, temp, nelq, nocs)
    net_Density_min = -netCharge(D_ef_min, ss) + net_Doping
    net_Density_max = -netCharge(D_ef_max, ss) + net_Doping
    
    # extend initial interval if necessary
    ei = 0
    while (net_Density_max < 0 and net_Density_min < 0) or (net_Density_max > 0 and net_Density_min > 0):
        if net_Density_max < 0 and net_Density_min < 0:
            ef_min = ef_min - (ef_max - ef_min)
            #print 'Interval adjusted. ef_min new = ', ef_min 
        if net_Density_max > 0 and net_Density_min > 0:
            ef_max = ef_max + (ef_max - ef_min)
            #print 'Interval adjusted. ef_max new = ', ef_max
        # calculate desities with new ef values
        D_ef_min = eDensity(WF_array, m_eff_array_q, EV, ef_min, temp, nelq, nocs)
        D_ef_max = eDensity(WF_array, m_eff_array_q, EV, ef_max, temp, nelq, nocs)
        net_Density_min = -netCharge(D_ef_min, ss)  + net_Doping
        net_Density_max = -netCharge(D_ef_max, ss)  + net_Doping
        ei += 1
        if ei > 100:
            intervalok = False
            break

    if intervalok == False:
        print 'could not find right interval!'
        return 0    
    
    ii = 0
    while abs(net_Density_min - net_Density_max) > net_Doping*1e-12:
        # calculate density for fermi level between ef_min and ef_max
        ef_middle = (ef_max + ef_min)/2
        D_ef_middle = eDensity(WF_array, m_eff_array_q, EV, ef_middle, temp, nelq, nocs)
        net_Density_middle = -netCharge(D_ef_middle, ss)  + net_Doping
        #print net_Density_middle
        if net_Density_middle < 0:
            ef_max = ef_middle
            D_ef_max = eDensity(WF_array, m_eff_array_q, EV, ef_max, temp, nelq, nocs)
            net_Density_max = -netCharge(D_ef_max, ss)  + net_Doping
        elif net_Density_middle > 0:
            ef_min = ef_middle
            D_ef_min = eDensity(WF_array, m_eff_array_q, EV, ef_min, temp, nelq, nocs)
            net_Density_min = -netCharge(D_ef_min, ss)  + net_Doping
        #print 'delta E = ', ef_max - ef_min
        #print 'ef_max = ', ef_max
        #print abs(net_Density_min - net_Density_max)
        if net_Density_min < 0 or net_Density_max > 0:
            print 'error'
            break
        ii += 1

    # check if converged
    D_ef_final = eDensity(WF_array, m_eff_array_q, EV, ef_max, temp, nelq, nocs)
    net_Density_final = -netCharge(D_ef_final, ss)  + net_Doping
    #print net_Density_final
    #print ii
    #print ef_max
    #print ef_min
    print 'new fermi level =', ef_max
    return ef_max
    
def V_ex(eDens_array, epsilon_array_q, m_eff_array_q):
    """ Exchange correlation term """
    
    a_star = 4.0*np.pi*const.epsilon_0*epsilon_array_q*(const.hbar**2)/(const.m_e*m_eff_array_q*const.e**2)
    r_s_inv = np.where(eDens_array>0., ((4.0/3*np.pi*a_star**3*eDens_array)**(1./3.)), 0.)
    Ry_star = (const.e**2)/(8.0*np.pi*const.epsilon_0*epsilon_array_q*a_star)
    V_ex = (-2.0*Ry_star)/(np.pi**(2.0/3)*(4.0/9)**(1.0/3))*r_s_inv + \
            (-2.0*Ry_star)/(np.pi**(2.0/3)*(4.0/9)**(1.0/3))*(0.7734/21)*np.log(1.0 + 21.0*r_s_inv)
    
    return V_ex

def plot_output(x, x_q, pot_tot_array_p, doping_n_array, eDens_array, nel, ef, time_ex, noit, target_error_p, 
                error_p, target_error_d, error_d, nocs, E, PSI, ss, gs, nomaxit, exchange_correlation_term, DEBUG, DEBUG_level, 
                fraction_in_dx_centers, fraction_of_free_charges_on_surface, surface_charge_array, temp, linear_solver, surface_charge_on):
    """ Output all relevant information and write into file """
    if DEBUG and (DEBUG_level==1 or DEBUG_level==2):
        plt.figure()
        plt.plot(x, pot_tot_array_p)
        for k in xrange(nocs):
            print "E[" + str(k) + "]=" + str(E[k])
            plt.plot(x_q, PSI[k]/np.max(abs(PSI[k])) + E[k])
        plt.title('Conduction band, wave functions, fermi level')    
        plt.show()
    
    # calculate sigma_z for file name
    s_z = sigma_z(x_q, eDens_array)
    if exchange_correlation_term:
        myfile = open('output_ex_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.txt', 'w')
    else: myfile = open('output_'  + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.txt', 'w')
    doping_sheet_density_cm = netCharge(doping_n_array, ss)*1e-4
    myfile.write('doping sheet density = ' + str('%.6e' % doping_sheet_density_cm) + ' cm^-2\n')  # convert from m^-2 to cm^-2 and write to output file
    electron_sheet_density_cm = -netCharge(eDens_array, ss)*1e-4
    myfile.write('electron sheet density = ' + str('%.6e' % electron_sheet_density_cm) + ' cm^-2\n') 
    surface_charge_sheet_density_cm = netCharge(surface_charge_array, ss)*1e-4
    myfile.write('surface charge sheet density = ' + str('%.6e' % surface_charge_sheet_density_cm) + ' cm^-2\n') 
    total_sheet_density_cm = (netCharge(doping_n_array, ss)-netCharge(eDens_array, ss)+netCharge(surface_charge_array, ss))*1e-4
    myfile.write('total sheet density = ' + str('%.6e' % total_sheet_density_cm) + ' cm^-2\n')
    myfile.write('fraction_in_dx_centers = ' + str(fraction_in_dx_centers) + '\n')
    myfile.write('fraction_of_free_charges_on_surface = ' + str(fraction_of_free_charges_on_surface) + '\n')
    myfile.write('eigenvalues = ' + str(E) + ' eV\n')
    myfile.write('final weights of wavefunctions = ' + str(wf_weights(PSI, E, nocs, ef, temp)) + '\n')
    myfile.write('converged at step ' + str(noit) + '\n')
    myfile.write('max number of iterations = ' + str(nomaxit) + '\n')
    myfile.write('target error poisson = ' + str(target_error_p) + ' V\n')
    myfile.write('error poisson = ' + str(error_p) + ' V\n')
    myfile.write('target error electron density = ' + str(target_error_d) + ' m^-3\n')
    myfile.write('error electron density = ' + str(error_d) + ' m^-3\n')
    myfile.write('number of elements = ' + str(nel) + '\n')
    myfile.write('number of considered states of schroedinger equation = ' + str(nocs) + '\n')
    myfile.write('exchange correlation = ' + str(exchange_correlation_term) + '\n')
    myfile.write('calculation time = ' + str(time_ex) + ' sec\n')
    myfile.write('fermi level = ' + str(ef) + ' eV\n')
    myfile.write('grid spacing = ' + str(gs) + ' nm\n')
    myfile.write('sigma_z = ' + str(s_z) + ' nm\n')
    myfile.write('temperature = ' + str(temp) + ' K\n')
    myfile.write('linear solver = ' + str(linear_solver) + '\n')
    myfile.write('surface charge = ' + str(surface_charge_on))
    myfile.close()
    
    # save file with electron density array in 1e24 m^-3 or equivalent 1e18 cm^-3
    if exchange_correlation_term:
        np.savetxt('eDens_array_ex_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) +  '_' + str(round(s_z,3)) + '.out', eDens_array*1e-24)
    else: np.savetxt('eDens_array_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) +  '_' + str(round(s_z,3)) + '.out', eDens_array*1e-24)
    if exchange_correlation_term:
        np.savetxt('pot_tot_array_p_ex_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.out', pot_tot_array_p)
    else: np.savetxt('pot_tot_array_p_' + str(nel) + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.out', pot_tot_array_p)
    if exchange_correlation_term:
        np.savetxt('Psi_ex_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.out', np.hstack([x_q,PSI.T]))
    else: np.savetxt('Psi_' + str(nel)  + '_' + str(fraction_in_dx_centers) + '_' + str(round(fraction_of_free_charges_on_surface,3)) + '_' + str(round(s_z,3)) + '.out', np.hstack([x_q,PSI.T]))
    print 'doping net charge = ', netCharge(doping_n_array, ss)
    print 'eDensity net charge = ', netCharge(eDens_array, ss)
    print 'total net charge = ', netCharge(doping_n_array, ss)-netCharge(eDens_array, ss)+netCharge(surface_charge_array, ss)
    print 'fraction_in_dx_centers = ', fraction_in_dx_centers
    print 'fraction_of_free_charges_on_surface = ', fraction_of_free_charges_on_surface
    print 'eigenvalues = ', E
    print 'final weights of wavefunctions = ', wf_weights(PSI, E, nocs, ef, temp)
    print 'converged at step', noit
    print 'target error_p =', target_error_p
    print 'error_p =', error_p
    print 'target error_d =', target_error_d
    print 'error_d =', error_d
    print 'number of elements =', nel
    print 'number of considered states of schroedinger equation =', nocs
    print 'exchange correlation =', exchange_correlation_term
    print 'sigma_z =', s_z
    print 'temperature = ', temp
    print 'linear solver = ', linear_solver
    print 'surface charge = ', surface_charge_on
    print 'calculation time =', time_ex

