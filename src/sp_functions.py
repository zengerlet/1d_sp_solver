import numpy as np
import scipy.integrate as integ
import scipy.constants as const
import matplotlib.pyplot as plt

def normalize(f, ss):
	'Normalize wave function'
	N = integ.simps(f*f, dx=ss)
	c = 1.0/np.sqrt(N)
	f_norm = c*f
	return f_norm

def netCharge(DOP_array, ss):
	'Integrate charge density array in z-direction'
	DOP_net_charge = integ.simps(DOP_array, dx=ss) 
	return DOP_net_charge

def eDensity(WF_array, m_eff_array_q, EV, ef, temp, nelq, nocs, oW=False):
	'Calculate eDensity from wave functions considering fermi-dirac distribution and density of states'
	density = np.zeros(nelq + 1)

	for i in xrange(nocs-1,-1,-1):
		exparg = -(EV[i] - ef)*const.e/const.k/temp
		# approximate for large arguments
		if exparg > 100.:
			expof = exparg
		else: expof = np.log(1 + np.exp(exparg))
		if oW:
			print 'Weight of WF', i+1, expof
		density = density + WF_array[i]**2*expof
	# multiply by constants and density of states 
	DOS = m_eff_array_q*const.m_e/const.pi/const.hbar**2
	density = density*temp*const.k*DOS
	return density
    
def findFermi(net_Doping, WF_array, EV, m_eff_array_q, temp, ss, nelq, nocs):
	'Find fermi level using charge neutrality'
	intervalok = True
	#start with ef_min = 0.0, ef_max = 1.0
	ef_min = -1.0
	ef_max = 1.0
	# check if given interval is ok
	D_ef_min = eDensity(WF_array, m_eff_array_q, EV, ef_min, temp, nelq, nocs)
	D_ef_max = eDensity(WF_array, m_eff_array_q, EV, ef_max, temp, nelq, nocs)
	net_Density_min = -netCharge(D_ef_min, ss) + net_Doping
	net_Density_max = -netCharge(D_ef_max, ss) + net_Doping
	
	# extend interval if necessary
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
	'Exchange correlation term'
	
	a_star = 4.0*np.pi*const.epsilon_0*epsilon_array_q*(const.hbar**2)/(const.m_e*m_eff_array_q*const.e**2)
	r_s_inv = np.where(eDens_array>0., ((4.0/3*np.pi*a_star**3*eDens_array)**(1./3.)), 0.)
	Ry_star = (const.e**2)/(8.0*np.pi*const.epsilon_0*epsilon_array_q*a_star)
	V_ex = (-2.0*Ry_star)/(np.pi**(2.0/3)*(4.0/9)**(1.0/3))*r_s_inv + \
			(-2.0*Ry_star)/(np.pi**(2.0/3)*(4.0/9)**(1.0/3))*(0.7734/21)*np.log(1.0 + 21.0*r_s_inv)
	
	'''
	plt.figure()
	plt.plot(V_ex/const.e)
	plt.title('exchange correlation potential')
	plt.show()
	#print 1/r_s_inv
	'''
	return V_ex

'''
def plot_output(x, x_q, pot_tot_array, doping_n_array, eDens_array, nel, ef, time_ex, noit, target_error, error, nocs, E, PSI, ss, nomaxit):
	'Output all relevant information and write into file'
	plt.figure()
	plt.plot(x_q, pot_tot_array)
	for k in xrange(nocs):
		print "E[" + str(k) + "]=" + str(E[k])
		plt.plot(x_q, PSI[k]/np.max(abs(PSI[k])) + E[k])
	plt.title('Conduction band, wave functions, fermi level')	
	plt.figure()
	plt.plot(E,'x')
	plt.title('Eigenvalues')
	plt.show()
	myfile = open('output_' + str(nel) + '.txt', 'w')
	myfile.write('doping net charge = ' + str(netCharge(doping_n_array, ss)) + '\n')
	myfile.write('eDensity net charge = ' + str(netCharge(eDens_array, ss)) + '\n') 
	myfile.write('eigenvalues = ' + str(E) + '\n')
	myfile.write('converged at step ' + str(noit) + '\n')
	myfile.write('max number of iterations = ' + str(nomaxit) + '\n')
	myfile.write('target error = ' + str(target_error) + '\n')
	myfile.write('error = ' + str(error) + '\n')
	myfile.write('number of elements = ' + str(nel) + '\n')
	myfile.write('number of considered states of schroedinger equation = ' + str(nocs) + '\n')
	myfile.write('calculation time = ' + str(time_ex) + 'sec\n')
	myfile.write('fermi level = ' + str(ef))
	myfile.close()
	np.savetxt('x_array_' + str(nel) + '.out',x)
	np.savetxt('x_q_array_' + str(nel) + '.out',x_q)
	np.savetxt('eDens_array_' + str(nel) + '.out', eDens_array)
	np.savetxt('pot_tot_array_' + str(nel) + '.out', pot_tot_array)
	print 'doping net charge = ', netCharge(doping_n_array, ss)
	print 'eDensity net charge = ', netCharge(eDens_array, ss)
	print 'converged at step', noit
	print 'target error =', target_error
	print 'error =', error
	print 'number of elements =', nel
	print 'number of considered states of schroedinger equation =', nocs
	print 'calculation time =', time_ex
'''

