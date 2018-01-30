# python 3 version
import matplotlib.pyplot as plt
import numpy as np
# Benjamin Klimko, PHYS 416 Spring 2018
#Code originally by Frank Toffoletto, edited by B. Klimko
# The program calculates the height of two objects numerically and then based on theory and plots all

# define interpolation functions
def intrpf(xi, x, y):
	yi = 0
	for q in range(0, len(y)):
		yi = yi + mult(xi, x, y, q)
	return yi

def mult(xi, x, y, q):
	num = 1.0
	denom = 1.0
	xnew = np.delete(x,q)
	for p in xnew:
		num = num*(xi-p)
		denom = denom*(x[q]-p)
	y_out = num*y[q]/denom
	return y_out

def balle_2(y1,speed, theta):  # compute two objects at once
	r1 = np.array([0.0, y1]);     # Initial vector position
	r2 = np.array([0.0, y1])
	v1 = np.array([speed*np.cos(theta*np.pi/180), 
	speed*np.sin(theta*np.pi/180)])  # Initial velocity
	v2 = np.array([speed*np.cos(theta*np.pi/180), 
	speed*np.sin(theta*np.pi/180)])  # Initial velocity
	r = np.copy(r1)
	v = np.copy(v1)  # Set initial position and velocity, best to copy to avoid overwrites
	r_2 = np.copy(r2)
	v_2 = np.copy(v2)
	time_mat = np.array([0.0])
	#* Set physical parameters (mass, Cd, etc.)
	Cd   = 0.5;    # Drag coefficient (dimensionless)
	big_area = 0.039;  # Cross-sectional area of projectile (m^2)
	small_area = 0.00182
	grav = 9.81;    # Gravitational acceleration (m/s^2)
	big_mass = 45.455;   # Mass of projectile (kg)  100 lbs
	small_mass = 0.455; # 1 lb

	rho = 1.2;    # Density of air (kg/m^3)
	air_const = -0.5*Cd*rho*big_area/big_mass;  # Air resistance constant
	air_const2 = -0.5*Cd*rho*small_area/small_mass
	#* Loop until ball hits ground or max steps completed
	tau = 0.01
	maxstep = 10000;   # Maximum number of steps
	for istep in range(0,maxstep):
		#* Record position (computed and theoretical) for plotting
		t = (istep)*tau;     # Current time
		time_mat = np.append(time_mat, t)
		if(istep ==0):
			xplot = np.array(r[0]);   # Record trajectory for plot
			yplot = np.array(r[1]);
			xplot_2 = np.array(r_2[0])
			yplot_2 = np.array(r_2[1])
		else:
			xplot = np.append(xplot,r[0]);   # Record trajectory for plot
			yplot = np.append(yplot,r[1]);
			xplot_2 = np.append(xplot_2,r_2[0]);   # Record trajectory for plot
			yplot_2 = np.append(yplot_2,r_2[1]);
			#* Calculate the acceleration of the ball
		accel = air_const*np.linalg.norm(v)*v;   # Air resistance
		accel_2 = air_const2*np.linalg.norm(v_2)*v_2
		accel[1] = accel[1]-grav;      # Gravity
		accel_2[1] = accel_2[1]-grav
		#* Calculate the new position and velocity using midpoint method
		v_old = v
		v_old2 = v_2 
		v = v + tau*accel;
		r = r + tau*(v+v_old)/2   
		v_2 = v_2 + tau*accel_2;
		r_2 = r_2 + tau*(v_2+v_old2)/2              
		
		#* If ball reaches ground (y<0), break out of the loop
		if( r[1] < 0 ):
			xplot = np.append(xplot,r[0])   # Record trajectory for plot
			yplot = np.append(yplot,r[1])
			xplot_2 = np.append(xplot_2,r_2[0])   # Record trajectory for plot
			yplot_2 = np.append(yplot_2,r_2[1])
			break                 # Break out of the for loop
	################### CODE FOR PART B
	b_heavy = (Cd*rho*big_area)/(2*big_mass)
	b_light = (Cd*rho*small_area)/(2*small_mass)
	theory_ylight = np.array([])
	theory_yheavy = np.array([])
	for element in time_mat:
		theory_ylight = np.append(theory_ylight, r2[1]- (1/b_light)*np.log(np.cosh(np.sqrt(b_light*grav)*element)))
		theory_yheavy = np.append(theory_yheavy, r1[1]- (1/b_heavy)*np.log(np.cosh(np.sqrt(b_heavy*grav)*element)))

	return time_mat, yplot, yplot_2, theory_yheavy, theory_ylight


# set values for initial height, speed, etc
y_i = 50
v_o = 0
theta_0 = 0


times, y_heavy, y_light, theo_h, theo_l = balle_2(y_i, v_o, theta_0)
# interpolate back to find when the 100 lb ball hit the ground (not below)
t_ground = intrpf(0, y_heavy[-3:], times[-3:])
light_height = intrpf(t_ground, times[-3:], y_light[-3:])

plt.figure()
plt.plot(times, y_heavy, label='100 lb ball')
plt.plot(times, y_light, label='1 lb ball')
plt.plot(times, theo_h, label='Theoretical 100 lb trajectory')
plt.plot(times, theo_l, label='Theoretical 1 lb trajectory')
plt.xlabel('Time (s)')
plt.ylabel('Height (m)')
plt.legend()
plt.title('Height over time for different ball weights')
plt.show()

print('When the heavier ball hits the ground the two are ', light_height, ' m apart')
if ((light_height)*100)/2.54 >= 2:
	print('Sadly for Galileo the small ball is not 2 inches away but ', ((y_light[-1]-y_heavy[-1])*100)/2.54, ' inches')
else:
	print('Galileo wins again, the wily bastard')