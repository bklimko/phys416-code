# python 3 version
import matplotlib.pyplot as plt
import numpy as np
# Program originally by Frank Toffoletto, modified by Benjamin Klimko
# for PHYS 416 Spring 2018

# define interpolation functions
def intrpf(xi, x, y):
	yi = 0
	for q in range(0, len(y)):
		yi = yi + mult(xi, x, y, q)
	return yi

def mult(xi, x, y, q):
	num = 1
	denom = 1
	xnew = np.delete(x,q)
	for p in xnew:
		num = num*(xi-p)
		denom = denom*(x[q]-p)
	y_out = num*y[q]/denom
	return y_out

#* Set initial position and velocity of the baseball
y1 = float(input("Enter initial height (meters): "))
r1 = np.array([0.0, y1])     # Initial vector position
t_mat = np.array([]) # time vector
kin = np.array([])  # kinetic energy vector
pot = np.array([])  # potential energy vector
tot = np.array([])  # total energy vector

speed = float(input("Enter initial speed (m/s): "))
theta = float(input("Enter initial angle (degrees): "))
v1 = np.array([speed*np.cos(theta*np.pi/180), 
speed*np.sin(theta*np.pi/180)])  # Initial velocity
r = np.copy(r1)
v = np.copy(v1)  # Set initial position and velocity, best to copy to avoid overwrites
#* Set physical parameters (mass, Cd, etc.)
Cd = 0.35   # Drag coefficient (dimensionless)
area = 4.3e-3  # Cross-sectional area of projectile (m^2)
grav = 9.81    # Gravitational acceleration (m/s^2)
mass = 0.145   # Mass of projectile (kg)
airFlag = float(input("Air resistance? (Yes:1, No:0):"))
if (airFlag == 0 ):
	rho = 0      # No air resistance
else:
	rho = 1.2    # Density of air (kg/m^3)
air_const = -0.5*Cd*rho*area/mass  # Air resistance constant
#* Loop until ball hits ground or max steps completed
tau = float(input("Enter timestep, tau (sec): "))  # (sec)
maxstep = 10000   # Maximum number of steps
for istep in range(0,maxstep):
	#* Record position (computed and theoretical) for plotting
	t = (istep)*tau     # Current time
	t_mat = np.append(t_mat, t)
	if(istep ==0):
		xplot = np.array(r[0])   # Record trajectory for plot
		yplot = np.array(r[1])
		xNoAir = np.array(r[0])
		yNoAir = np.array(r[1])
	else:
		xplot = np.append(xplot,r[0])   # Record trajectory for plot
		yplot = np.append(yplot,r[1])
		xNoAir = np.append(xNoAir,r1[0] + v1[0]*t)   # Record trajectory for plot
		yNoAir = np.append(yNoAir,r1[1] + v1[1]*t - 0.5*grav*t**2)
	#* Calculate the acceleration of the ball
	accel = air_const*np.linalg.norm(v)*v  # Air resistance
	accel[1] = accel[1]-grav      # Gravity
	# calculate the kinetic, potential, and total energies of the system
	kin = np.append(kin, .5 * mass * np.linalg.norm(v) ** 2)
	pot = np.append(pot, mass * grav * r[1])
	tot = np.append(tot, kin[-1]+pot[-1])
	#* Calculate the new position and velocity using Euler method
	r = r + tau*v                 # Euler step
	v = v + tau*accel

	#* If ball reaches ground (y<0), break out of the loop
	if( r[1] < 0 ):
		xplot = np.append(xplot,r[0])   # Record trajectory for plot
		yplot = np.append(yplot,r[1])
		break                  # Break out of the for loop
# Print maximum range and time of flight
x = intrpf(0, [yplot[-3], yplot[-2], yplot[-1]], [xplot[-3], xplot[-2], xplot[-1]])
t_r = intrpf(0, [yplot[-3], yplot[-2], yplot[-1]],[t_mat[-3], t_mat[-2], t_mat[-1]])
print('Maximum range is ',x,' meters')
print('Time of flight is ',t_r,' seconds')
# Graph the trajectory of the baseball
plt.figure(0);   plt.clf(); # Clear figure window and bring it forward
# Mark the location of the ground by a straight line
xground = np.array([0, np.max(xNoAir)]);  yground = np.array([0, 0]);
# Plot the computed trajectory and parabolic, no-air curve
plt.plot(xplot,yplot,'.')
plt.plot(xNoAir,yNoAir,'--')
plt.plot(xground,yground,'-')
plt.legend(['Euler method','Theory (No air)'])
plt.xlabel('Range (m)');  plt.ylabel('Height (m)')
plt.title('Projectile motion')
plt.savefig('projectile.png', bbox_inches='tight')

# plot the various energies
plt.figure(1)
plt.plot(t_mat, kin)
plt.xlabel('Time (s)')
plt.ylabel('Kinetic Energy (J)')
plt.title('Kinetic energy of system')
plt.savefig('kineticE.png', bbox_inches='tight')

plt.figure(2)
plt.plot(t_mat, pot)
plt.xlabel('Time (s)')
plt.ylabel('Potential Energy (J)')
plt.title('Potential energy of system')
plt.savefig('potentialE.png', bbox_inches='tight')

plt.figure(3)
plt.plot(t_mat, tot)
plt.xlabel('Time (s)')
plt.ylabel('Total Energy (J)')
plt.title('Total energy of system')
plt.ylim(180, 190)
plt.savefig('totalE.png', bbox_inches='tight')
#axis equal; shg; # reset the aspect ratio, bring the plot to the front
plt.grid('on')
plt.show()