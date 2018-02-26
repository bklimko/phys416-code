import numpy as np
import matplotlib.pyplot as plt
def rk4(x,t,tau,derivsRK, scale):
	#%  Runge-Kutta integrator (4th order)
	#% Input arguments -
	#%   x = current value of dependent variable
	#%   t = independent variable (usually time)
	#%   tau = step size (usually timestep)
	#%   derivsRK = right hand side of the ODE; derivsRK is the
	#%             name of the function which returns dx/dt
	#%             Calling format derivsRK(x,t).
	#% Output arguments -
	#%   xout = new value of x after a step of size tau
	half_tau = 0.5*tau
	F1 = derivsRK(x,t, scale)
	t_half = t + half_tau
	xtemp = x + half_tau*F1
	F2 = derivsRK(xtemp,t_half, scale)
	xtemp = x + half_tau*F2
	F3 = derivsRK(xtemp,t_half, scale)
	t_full = t + tau
	xtemp = x + tau*F3
	F4 = derivsRK(xtemp,t_full, scale)
	xout = x + tau/6.*(F1 + F4 + 2.*(F2+F3))
	return xout
def rka(x,t,tau,err,derivsRK, scale):
	#% Adaptive Runge-Kutta routine
	#% Inputs
	#%   x          Current value of the dependent variable
	#%   t          Independent variable (usually time)
	#%   tau        Step size (usually time step)
	#%   err        Desired fractional local truncation error
	#%   derivsRK   Right hand side of the ODE; derivsRK is the
	#%              name of the function which returns dx/dt
	#%              Calling format derivsRK(x,t).
	#% Outputs
	#%   xSmall     New value of the dependent variable
	#%   t          New value of the independent variable
	#%   tau        Suggested step size for next call to rka
	#%* Set initial variables
	tSave = t;  xSave = x    # Save initial values
	safe1 = .9;  safe2 = 4.  # Safety factors
	eps = np.spacing(1) # smallest value
	#%* Loop over maximum number of attempts to satisfy error bound
	maxTry = 100
	for iTry in range(1,maxTry):
		#%* Take the two small time steps
		half_tau = 0.5 * tau
		xTemp = rk4(xSave,tSave,half_tau,derivsRK, scale)
		t = tSave + half_tau
		xSmall = rk4(xTemp,t,half_tau,derivsRK, scale)
		#%* Take the single big time step
		t = tSave + tau
		xBig = rk4(xSave,tSave,tau,derivsRK, scale)
		#%* Compute the estimated truncation error
		scale = err * (np.abs(xSmall) + np.abs(xBig))/2.
		xDiff = xSmall - xBig
		errorRatio = np.max( [np.abs(xDiff)/(scale + eps)] )
		#print safe1,tau,errorRatio
		#%* Estimate news tau value (including safety factors)
		tau_old = tau
		tau = safe1*tau_old*errorRatio**(-0.20)
		tau = np.max([tau,tau_old/safe2])
		tau = np.min([tau,safe2*tau_old])
		#%* If error is acceptable, return computed values
		if errorRatio < 1 :
			# xSmall = xSmall #% +  (xDiff)/15
			#   xSmall = (16.*xSmall - xBig)/15. # correction
			return xSmall, t, tau
	#%* Issue error message if error bound never satisfied
	print ('ERROR: Adaptive Runge-Kutta routine failed')
	return
def gravrk(s,t, scale):
	#%  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
	#%  Inputs
	#%    s      State vector [r(1) r(2) v(1) v(2)]
	#%    t      Time (not used)
	#     r_other     position of other body
	#     scale       Scaling weight of other body relative to the sun
	#%  Output
	#%    deriv  Derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt] (expanded out for two bodies)
	
	r = np.array([s[0], s[1]])  # Unravel the vector s into position and velocity
	v = np.array([s[2] ,s[3]])
	r_other = np.array([s[4], s[5]])
	v_other = np.array([s[6], s[7]])

	# compute accelerations for both bodies
	accel =  (-GM*r/(np.linalg.norm(r)**3)) -scale[1]*(GM*(r - r_other))/(np.linalg.norm(r - r_other)**3) # Gravitational acceleration
	accel_other =  (-GM*r_other/(np.linalg.norm(r_other)**3)) -scale[0]*(GM*(r_other - r))/(np.linalg.norm(r_other - r)**3)
	 
	
	#%* Return derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
	derivs = np.array([v[0], v[1], accel[0], accel[1], v_other[0], v_other[1], accel_other[0], accel_other[1]])
	return derivs

# orbit - Program to compute the orbits of a Sun-Earth-Jupiter system with various masses of Jupiter
# Benjamin Klimko, PHYS 416 Spring 2018

# Set initial position and velocity of the planets.
r_e = np.array([1., 0.]) # earth and jupiter position vectors  
r_j = np.array([5.2, 0.])
v_e = np.array([0., 2*np.pi]) # earth and jupiter velocities
v_j = np.array([0., 2*np.pi/np.sqrt(5.2)])

state = np.array([ r_e[0], r_e[1], v_e[0], v_e[1], r_j[0], r_j[1], v_j[0], v_j[1] ])   # Used by R-K routines

#Set physical parameters (mass, G*M)
GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)

masses = np.array([6.0e24, 1000*1.9e27])/2e30
adaptErr = 1.e-4 # Error parameter used by adaptive Runge-Kutta
time = 0.0


#%* Loop over desired number of steps using adaptive Runge-Kutta
nStep = int(input("Enter number of steps: "))
tau =  float(input("Enter time step (yr): "))

for istep in range(0,nStep):
	#%* Record position and energy for plotting.
	# Initially set the arrays for the first step
	if istep == 0:
		rplot_e = np.linalg.norm(r_e)
		thplot_e = np.arctan2(r_e[1],r_e[0])
		tplot = time
		rplot_j = np.linalg.norm(r_j)
		thplot_j = np.arctan2(r_j[1],r_j[0])
		
	else:
		rplot_e = np.append(rplot_e,np.linalg.norm(r_e))           #Record position for earth
		thplot_e = np.append(thplot_e,np.arctan2(r_e[1],r_e[0]))
		tplot = np.append(tplot,time)
		rplot_j = np.append(rplot_j,np.linalg.norm(r_j))           #Record position for Jupiter
		thplot_j = np.append(thplot_j,np.arctan2(r_j[1],r_j[0]))
		
	#%* Calculate new position and velocity using adaptive Runge-Kutta
	[state, time, tau] = rka(state,time,tau,adaptErr,gravrk, masses)
	
	# separate out each body's position and velocity from the state vector
	r_e = np.array([state[0], state[1]] )
	v_e = np.array([state[2], state[3]])
	r_j =  np.array([state[4] ,state[5]])
	v_j = np.array([state[6], state[7]])

	
#%* Graph the trajectory of the comet.

xplot_e = rplot_e * np.cos(thplot_e)
yplot_e = rplot_e * np.sin(thplot_e)

xplot_j = rplot_j * np.cos(thplot_j)
yplot_j = rplot_j * np.sin(thplot_j)

plt.figure(1); plt.clf()  #Clear figure 1 window and bring forward

plt.plot(xplot_e,yplot_e,label='Earth Trajectory')  # Use plot for graphing 
plt.plot(xplot_j, yplot_j, label='Jupiter Trajectory')
plt.legend()
plt.title('Sun-Earth-Jupiter system')
plt.xlabel('Distance (AU)')
plt.grid(True)
plt.show()
