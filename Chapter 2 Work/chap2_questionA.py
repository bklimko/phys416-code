import matplotlib.pyplot as plt
import numpy as np
# Benjamin Klimko, PHYS 416 Spring 2018
#Code originally by Frank Toffoletto, edited by above
#* Set initial position and velocity of the baseball
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
	
	area = 0.002827;  # Cross-sectional area of projectile (m^2)
	grav = 9.81;    # Gravitational acceleration (m/s^2)
	mass = 0.04593;   # Mass of projectile (kg)

	rho = 1.2;    # Density of air (kg/m^3)
	Sw = [0, 0, 0.25]
	
	#* Loop until ball hits ground or max steps completed
	tau = 0.01
	maxstep = 10000;   # Maximum number of steps
	for istep in range(0,maxstep):
		# Drag coefficient (dimensionless)
		if np.linalg.norm(v) <= 14:
			Cd = 0.5
		else:
			Cd = 7.0/speed
		air_const = -0.5*Cd*rho*area/mass;  # Air resistance constant
		magnus = np.cross(Sw, [v[0], v[1], 0])
		
		#* Record position (computed and theoretical) for plotting
		t = (istep)*tau;     # Current time
		if(istep ==0):
			xplot = np.array(r[0]);   # Record trajectory for plot
			yplot = np.array(r[1]);
			
		else:
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
		
			#* Calculate the acceleration of the ball
		accel = air_const*np.linalg.norm(v)*v + magnus[:-1]   # Air resistance
	
		accel[1] = accel[1]-grav       # Gravity 
		#* Calculate the new position and velocity using midpoint method
		v_old = v 
		v = v + tau*accel;
		r = r + tau*(v+v_old)/2;                 
		
		#* If ball reaches ground (y<0), break out of the loop
		if( r[1] < 0 ):
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			break                 # Break out of the for loop
	drange = intrpf(0, yplot[-3:], xplot[-3:])
	return drange, xplot, yplot

def balle_smooth(y1,speed, theta):
	r1 = np.array([0.0, y1]);     # Initial vector position

	v1 = np.array([speed*np.cos(theta*np.pi/180), 
	speed*np.sin(theta*np.pi/180)])  # Initial velocity
	r = np.copy(r1)
	v = np.copy(v1)  # Set initial position and velocity, best to copy to avoid overwrites
	#* Set physical parameters (mass, Cd, etc.)
	# # Drag coefficient (dimensionless)
	Cd = 0.5
	area = 0.002827;  # Cross-sectional area of projectile (m^2)
	grav = 9.81;    # Gravitational acceleration (m/s^2)
	mass = 0.04593;   # Mass of projectile (kg)

	rho = 1.2;    # Density of air (kg/m^3)
	air_const = -0.5*Cd*rho*area/mass;  # Air resistance constant
	Sw = [0, 0, 0.25]
	#* Loop until ball hits ground or max steps completed
	tau = 0.01
	maxstep = 10000;   # Maximum number of steps
	for istep in range(0,maxstep):
		
		magnus = np.cross(Sw, [v[0], v[1], 0]) #magnus force
		#* Record position (computed and theoretical) for plotting
		t = (istep)*tau;     # Current time
		if(istep ==0):
			xplot = np.array(r[0]);   # Record trajectory for plot
			yplot = np.array(r[1]);
			
		else:
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			
			#* Calculate the acceleration of the ball
		accel = air_const*np.linalg.norm(v)*v+ magnus[:-1]   # Air resistance
		
		accel[1] = accel[1]-grav      # Gravity plus Magnus force
		#* Calculate the new position and velocity using midpoint method
		v_old = v 
		v = v + tau*accel;
		r = r + tau*(v+v_old)/2;                 
		
		#* If ball reaches ground (y<0), break out of the loop
		if( r[1] < 0 ):
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			break                 # Break out of the for loop
	srange = intrpf(0, yplot[-3:], xplot[-3:])
	return srange, xplot, yplot


# create vector of angles to test
angles = np.array([1, 3, 5, 7, 10])
angles_s = np.array([1,5,10, 15, 25])
# set values for initial height, speed, etc
y_i = 0
v_o = 70
ranges = np.array([])
ranges_smooth = np.array([])


#plot trajectories for both smooth and dimpled golf balls
plt.figure(1)
for angle in angles:
 	[range_ball, xvals, yvals] = balle(y_i, v_o, angle)
 	plt.plot(xvals, yvals, label=str(angle) + ' degrees')
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Trajectories of dimpled golf ball for various launch angles')
plt.legend()
plt.show()

plt.figure(2)
for angle in angles_s:
 	[range_ball, xvals, yvals] = balle_smooth(y_i, v_o, angle)
 	plt.plot(xvals, yvals, label=str(angle) + ' degrees')
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Trajectories of smooth golf ball for various launch angles')
plt.legend()
plt.show()

# find the optimal angle for maximum range (for both balls)
angles2 = np.linspace(0, 90)

# loop through every angle to find the range at that point for both balls
for angle2 in angles2:
	[r_d, x, y] = balle(y_i, v_o, angle2)
	[r_s, x1, y1] = balle_smooth(y_i, v_o, angle2)
	ranges = np.append(ranges, r_d)
	ranges_smooth = np.append(ranges_smooth, r_s)

max_range_dimple = np.amax(ranges)
max_angle_dimple = angles2[np.argmax(ranges)]
print('The maximum range is ', max_range_dimple, ' m at an angle of ', max_angle_dimple, ' degrees for a dimpled golf ball')
max_range_smooth = np.amax(ranges_smooth)
max_angle_smooth = angles2[np.argmax(ranges_smooth)]
print('The maximum range is ', max_range_smooth, ' m at an angle of ', max_angle_smooth, ' degrees for a smooth golf ball')