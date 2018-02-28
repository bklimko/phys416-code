import numpy as np 
import matplotlib.pyplot as plt 

# Program for problem 22 to calculate the trajectory of predator and prey populations
# Benjamin Klimko, PHYS 416, Spring 2018


def rk4(x,t,tau,derivsRK):
	##  Runge-Kutta integrator (4th order)
	## Input arguments -
	##   x = current value of dependent variable
	##   t = independent variable (usually time)
	##   tau = step size (usually timestep)
	##   derivsRK = right hand side of the ODE; derivsRK is the
	##             name of the function which returns dx/dt
	##             Calling format derivsRK(x,t).
	## Output arguments -
	##   xout = new value of x after a step of size tau
	half_tau = 0.5*tau
	F1 = derivsRK(x,t)
	t_half = t + half_tau
	xtemp = x + half_tau*F1
	F2 = derivsRK(xtemp,t_half)
	xtemp = x + half_tau*F2
	F3 = derivsRK(xtemp,t_half)
	t_full = t + tau
	xtemp = x + tau*F3
	F4 = derivsRK(xtemp,t_full)
	xout = x + tau/6.*(F1 + F4 + 2.*(F2+F3))
	return xout

def rka(x,t,tau,err,derivsRK):
	## Adaptive Runge-Kutta routine
	## Inputs
	##   x          Current value of the dependent variable
	##   t          Independent variable (usually time)
	##   tau        Step size (usually time step)
	##   err        Desired fractional local truncation error
	##   derivsRK   Right hand side of the ODE; derivsRK is the
	##              name of the function which returns dx/dt
	##              Calling format derivsRK(x,t).
	## Outputs
	##   xSmall     New value of the dependent variable
	##   t          New value of the independent variable
	##   tau        Suggested step size for next call to rka
	##* Set initial variables
	tSave = t;  xSave = x    # Save initial values
	safe1 = .9;  safe2 = 4.  # Safety factors
	eps = np.spacing(1) # smallest value
	##* Loop over maximum number of attempts to satisfy error bound
	maxTry = 100
	for iTry in range(1,maxTry):
		##* Take the two small time steps
		half_tau = 0.5 * tau
		xTemp = rk4(xSave,tSave,half_tau,derivsRK)
		t = tSave + half_tau
		xSmall = rk4(xTemp,t,half_tau,derivsRK)
		##* Take the single big time step
		t = tSave + tau
		xBig = rk4(xSave,tSave,tau,derivsRK)
		##* Compute the estimated truncation error
		scale = err * (np.abs(xSmall) + np.abs(xBig))/2.
		xDiff = xSmall - xBig
		errorRatio = np.max( [np.abs(xDiff)/(scale + eps)] )
		#print safe1,tau,errorRatio
		##* Estimate news tau value (including safety factors)
		tau_old = tau
		tau = safe1*tau_old*errorRatio**(-0.20)
		tau = np.max([tau,tau_old/safe2])
		tau = np.min([tau,safe2*tau_old])
		##* If error is acceptable, return computed values
		if errorRatio < 1 :
			xSmall = xSmall ## +  (xDiff)/15
			return xSmall, t, tau
	##* Issue error message if error bound never satisfied
	print ('ERROR: Adaptive Runge-Kutta routine failed')
	return

def lotkavolrk(s,t):
	#  Returns right-hand side of Lotka-Volterra model ODEs
	#  Inputs
	#    s      State vector [r f]
	#    t      Time (not used)
	#  Output
	#    deriv  Derivatives [dr/dt df/dt]
	#* For clarity, unravel input vectors
	r = s[0]; f = s[1];
	#* Return the derivatives [dr/dt df/dt]
	deriv = np.zeros(2)
	deriv[0] = a*(1-(r/b))*r - c*r*f 
	deriv[1] = -a*f + c*r*f 
	
	return deriv



# initialize constants and initial values
a = 10
b = 1e6
c = 0.1

r_0 = 100
f_0 = 100

state = np.array([r_0, f_0]) # state vector: r(t) f(t)

tau = 1       # Initial guess for the timestep
err = 1.e-3   # Error tolerance
time = 0

nstep = int(input('Enter number of steps: '))

# loop through each step
for istep in range(0, nstep):
	if istep == 0: # initialize data vectors
		rplot = np.array([state[0]])
		fplot = np.array([state[1]])
		tplot = np.array([time])

	else: # append to data vectors
		rplot = np.append(rplot, state[0])
		fplot = np.append(fplot, state[1])
		tplot = np.append(tplot, time)

	[state, time, tau] = rka(state, time, tau, err, lotkavolrk) # use adaptive Runge-Kutta to update the equations


# plot r(t) vs f(t)
plt.figure()
plt.plot(fplot, rplot, '.')
plt.xlabel('Fox Population')
plt.ylabel('Rabbit Population')
plt.show()