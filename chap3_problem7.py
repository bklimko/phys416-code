import numpy as np 
import matplotlib.pyplot as plt 

#Program to calculate trajectory of a charged particle in an electric and magnetic field
#Benjamin Klimko, PHYS 416 Spring 2018

#function to compute gyroradius
def gyro_r(B, mass, charge, vel):
	return (mass * np.linalg.norm(vel)) / (charge * np.linalg.norm(B))

# function to compute gyroperiod
def gyro_p(B, charge, mass):
	return (2* np.pi * mass) / (charge * np.linalg.norm(B))

# set normalized quantities
mass = 1
q = 1 # charge
# initial conditions
r = np.array([0., 0.])
v = np.array([0., 1.,0.])

# get user inputs for time step and total simulation time
tau = float(input('Enter time step (s): '))
totaltime = int(input('Enter total simulation time (s): '))

# set electric and magnetic field vectors
E = np.array([0.,0.,0.])
# B = np.array([0.,0.,1.])

# est_gr = gyro_r(B, mass, q, v)
# est_gp = gyro_p(B, q, mass)
# print('Estimated gyroradius: ', est_gr)
# print('Estimated gyroperiod: ', est_gp)

# calculate E x B drift
# v_eb = np.cross(E, B) / np.linalg.norm(B)**2


#determine how many time steps to take
# nSteps = int(np.ceil(est_gp / tau)) + 50  # for zero E field and constant B field
nSteps = int(np.ceil(totaltime/tau))

# for part A and B with constant B field
for step in range(0, nSteps):
	
	# calculate B and del-B drift
	B = np.array([0.,0., 1+(0.1*r[0])])
	v_delb = ((mass*np.linalg.norm(v)**2) / (2 * q * np.linalg.norm(B))) * (np.cross(B, np.gradient(B)) / np.linalg.norm(B)**2)

	# calculate acceleration
	accel = (q * (E + np.cross(v, B))) / mass

	if step == 0: # set up vectors of angle and radius data for plotting
		rplot = np.linalg.norm(r)
		thplot = np.arctan2(r[1],r[0])
		xplot = r[0]
		yplot = r[1]
		xplot_mf = r[0]
		yplot_mf = r[1]
		
		
	else: # append new data to data vectors
		rplot = np.append(rplot, np.linalg.norm(r))
		thplot = np.append(thplot, np.arctan2(r[1],r[0]))
		xplot = np.append(xplot, r[0])
		yplot = np.append(yplot, r[1])
		xplot_mf = np.append(xplot_mf, r[0]- (tau*v_delb[0]))
		yplot_mf = np.append(yplot_mf, r[1] - (tau*v_delb[1]))
		

	# Euler-Cromer step
	v = v + tau*accel
	r = r + tau*v[:-1]


# plotting in the x-y plane for rest frame and moving frame

# for E x B motion
# xplot_mf = xplot_mf - (tau*np.arange(len(xplot_mf))*v_eb[0])
# yplot_mf = yplot_mf - (tau*np.arange(len(yplot_mf))*v_eb[1])



plt.figure()
plt.plot(xplot, yplot)
plt.title('rest frame')
plt.grid(True)

plt.figure(2)
plt.plot(xplot_mf, yplot_mf)
plt.title('moving frame')
plt.grid(True)
plt.show()
print ('Actual gyroperiod: ', step*tau, 's')

