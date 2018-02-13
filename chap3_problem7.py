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
r_moveframe = np.copy(r) # create a copy of the radius to plot the trajectory in a moving reference frame

# get user inputs for time step and total simulation time
tau = float(input('Enter time step (s): '))
totaltime = int(input('Enter total simulation time (s): '))

# set electric and magnetic field vectors
E = np.array([1.,0.,0.])
B = np.array([0.,0.,1.])

est_gr = gyro_r(B, mass, q, v)
est_gp = gyro_p(B, q, mass)
print('Estimated gyroradius: ', est_gr)
print('Estimated gyroperiod: ', est_gp)

# calculate E x B drift as well as del-B drift
v_eb = np.cross(E, B) / np.linalg.norm(B)**2
v_delb = ((mass*np.linalg.norm(v)**2) / (2 * q * np.linalg.norm(B))) * (np.cross(B, np.gradient(B)) / np.linalg.norm(B)**2)

#determine how many time steps to take
# nSteps = int(np.ceil(est_gp / tau)) + 50
nSteps = int(np.ceil(totaltime/tau))

# for part A and B with constant B field
for step in range(0, nSteps):
	# calculate acceleration
	accel = (q * (E + np.cross(v, B))) / mass

	if step == 0:
		rplot = np.linalg.norm(r)
		thplot = np.arctan2(r[1],r[0])
		moveframe = np.linalg.norm(r_moveframe)
	else:
		rplot = np.append(rplot, np.linalg.norm(r))
		thplot = np.append(thplot, np.arctan2(r[1],r[0]))
		moveframe = np.append(moveframe, np.linalg.norm(r_moveframe))

	# Euler-Cromer step
	v = v + tau*accel
	r = r + tau*v[:-1]
	r_moveframe = r + tau*(v[:-1]-v_eb[:-1])


plt.figure()
plt.polar(thplot, rplot)
plt.title('stationary frame')
plt.grid(True)

plt.figure(2)
plt.polar(thplot, moveframe)
plt.title('moving frame')
plt.show()
print ('Actual gyroperiod: ', step*tau, 's')

