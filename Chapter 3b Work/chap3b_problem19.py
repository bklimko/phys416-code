import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  

# Program for problem 19 to calculate the trajectory of predator and prey populations
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

def dblpendrk(s,t):
	#  Returns right-hand side of Lotka-Volterra model ODEs
	#  Inputs
	#    s      State vector [theta1, theta2, dtheta1, dtheta2]
	#    t      Time (not used)
	#  Output
	#    deriv  Derivatives [dtheta1/dt, dtheta2/dt, d^2theta1/dt^2 d^2theta2/dt^2]
	
	alpha = (m1+m2)*L1
	beta = m2*L2*np.cos(s[0]-s[1])
	a = -(m1+m2)*g*np.sin(s[0]) - m2*L2*s[3]**2 *np.sin(s[0]-s[1])
	delta = m2*L1*np.cos(s[0]-s[1])
	ep = m2*L2
	b = -m2*g*np.sin(s[1]) + m2*L1*s[2]**2 * np.sin(s[0]-s[1])
	left = np.array([[alpha, beta], [delta, ep]])
	right = np.array([a, b])
	#* Return the derivatives [d^2theta1/dt^2 d^2theta2/dt^2]
	# dbldot = np.linalg.inv(left).dot(right)
	dbldot = np.linalg.solve(left, right)
	deriv = np.array([s[2], s[3], dbldot[0], dbldot[1]])
	
	return deriv

# set initial conditions and constants
g = 9.81 # gravity
m1 = 1
m2 = 1
L1 = 0.1
L2 = 0.1

#initial angles
theta1 = np.pi/4
theta2 = np.pi/2
th1dot = np.pi
th2dot = np.pi/12

# state vector: theta1, theta2, theta1_dot, theta2_dot
state = np.array([theta1, theta2, th1dot, th2dot]) 

tau = 1       # Initial guess for the timestep
err = 1.e-5  # Error tolerance
time = 0

nstep = 500

# loop through each step
for istep in range(0, nstep):
	if istep == 0: # initialize data vectors
		th1_plot = np.array([state[0]])
		th2_plot = np.array([state[1]])
		tplot = np.array([time])
		potential = np.array([-m1*L1*g*np.cos(state[0]) - m2*g*(L1*np.cos(state[0])+L2*np.cos(state[1]))])
		kinetic = np.array([(.5*(m1+m2)*L1**2*state[2]**2) + (.5*m2*L2**2*state[3]**2) + (m2*L1*L2*state[2]*state[3]*np.cos(state[0]-state[1]))])

	else: # append to data vectors
		th1_plot = np.append(th1_plot, state[0])
		th2_plot = np.append(th2_plot, state[1])
		tplot = np.append(tplot, time)
		potential = np.append(potential, -m1*L1*g*np.cos(state[0]) - m2*g*(L1*np.cos(state[0])+L2*np.cos(state[1])))
		kinetic = np.append(kinetic, (.5*(m1+m2)*L1**2*state[2]**2) + (.5*m2*L2**2*state[3]**2) +\
		 (m2*L1*L2*state[2]*state[3]*np.cos(state[0]-state[1])))

	[state, time, tau] = rka(state, time, tau, err, dblpendrk) # use adaptive Runge-Kutta to update the equations

# compute values for plotting

x1 = L1*np.sin(th1_plot)
y1 = -L1*np.cos(th1_plot)
x2 = x1 + L2*np.sin(th2_plot)
y2 = y1 - L2*np.cos(th2_plot)


#animate
plt.figure()
ax = plt.gca()
plt.title('Double Pendulum')

for i in range(0, len(x1)):
	plt.plot(x1[i], y1[i], 'ro')
	plt.plot(x2[i], y2[i], 'ko')
	ax.set_xlim(-.3, .3)
	ax.set_ylim(-.3, .3)
	l1 = Line2D([x1[i], x2[i]], [y1[i], y2[i]])
	l2 = Line2D([0, x1[i]], [0, y1[i]])
	ax.add_line(l1)
	ax.add_line(l2)
	plt.draw()
	plt.pause(0.01)
	# l1.remove()
	# l2.remove()
	if(i==len(x1)-1):
		#save plot out
		plt.savefig('c3bq19_3.png', bbox_inches='tight')
		# temp=input('Hit return to end:')
	else:
		# plt.clf()
		l1.remove()
		l2.remove()


totalE = potential + kinetic
plt.figure(2)
plt.plot(tplot, potential, label='potential energy')
plt.plot(tplot, kinetic, label='kinetic energy')
plt.plot(tplot, totalE, label='total energy')
plt.title('Double Pendulum Energy Diagram')
plt.legend()
plt.show()
# plt.savefig('c3bq19_1e.png', bbox_inches='tight')

temp=input('Hit return to end:')
