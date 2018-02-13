# python 3 version
import matplotlib.pyplot as plt
import numpy as np
# Benjamin Klimko, PHYS 416 Spring 2018
#Code originally by Frank Toffoletto, edited by above

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

def balle(y1,speed, theta):
	r1 = np.array([0.0, y1]);     # Initial vector position

	v1 = np.array([speed*np.cos(theta*np.pi/180), 
	speed*np.sin(theta*np.pi/180)])  # Initial velocity
	r = np.copy(r1)
	v = np.copy(v1)  # Set initial position and velocity, best to copy to avoid overwrites
	#* Set physical parameters (mass, Cd, etc.)
	Cd   = 0.35;    # Drag coefficient (dimensionless)
	area = 4.3e-3;  # Cross-sectional area of projectile (m^2)
	grav = 9.81;    # Gravitational acceleration (m/s^2)
	mass = 0.145;   # Mass of projectile (kg)

	rho = 1.2;    # Density of air (kg/m^3)
	air_const = -0.5*Cd*rho*area/mass;  # Air resistance constant
	#* Loop until ball hits ground or max steps completed
	tau = 0.01
	maxstep = 10000;   # Maximum number of steps
	for istep in range(0,maxstep):
		#* Record position (computed and theoretical) for plotting
		t = (istep)*tau;     # Current time
		if(istep ==0):
			xplot = np.array(r[0]);   # Record trajectory for plot
			yplot = np.array(r[1]);
			# xNoAir = np.array(r[0]);
			# yNoAir = np.array(r[1]);
		else:
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			# xNoAir = np.append(xNoAir,r1[0] + v1[0]*t);   # Record trajectory for plot
			# yNoAir = np.append(yNoAir,r1[1] + v1[1]*t - 0.5*grav*t**2);
			#* Calculate the acceleration of the ball
		accel = air_const*np.linalg.norm(v)*v;   # Air resistance
		accel[1] = accel[1]-grav;      # Gravity
		#* Calculate the new position and velocity using midpoint method
		v_old = v 
		v = v + tau*accel;
		r = r + tau*(v+v_old)/2;                 
		
		#* If ball reaches ground (y<0), break out of the loop
		if( r[1] < 0 ):
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			break                 # Break out of the for loop
	out = intrpf(0, yplot[-3:], xplot[-3:])
	return out 

# create vector of angles to test
angles = np.linspace(10, 50, 1000)
# set values for initial height, speed, etc
y_i = 1
v_o = 50
ranges = np.array([])


# loop through every angle to find the range at that point
for angle in angles:
	ranges = np.append(ranges, balle(y_i, v_o, angle))

max_range = np.amax(ranges)
max_angle = angles[np.argmax(ranges)]
print('The maximum range is ', max_range, ' m at an angle of ', max_angle, ' degrees')

plt.figure()
plt.plot(angles, ranges)
plt.xlabel(r'$Angle (\theta)$')
plt.ylabel('Range (m)')
plt.title('Range as a function of initial angle')
plt.show()