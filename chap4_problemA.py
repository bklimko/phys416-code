import matplotlib.pyplot as plt
import numpy as np
def rk4(x,t,tau,derivsRK):
	#%  Runge-Kutta integrator (4th order)
	#% Input arguments -
	#%   x = current value of dependent variable
	#%   t = independent variable (usually time)
	#%   tau = step size (usually timestep)
	#%   derivsRK = right hand side of the ODE; derivsRK is the
	#%             name of the function which returns dx/dt
	#% Output arguments -
	#%   xout = new value of x after a step of size tau
	half_tau = 0.5*tau;
	F1 = derivsRK(x,t);
	t_half = t + half_tau;
	xtemp = x + half_tau*F1;
	F2 = derivsRK(xtemp,t_half);
	xtemp = x + half_tau*F2;
	F3 = derivsRK(xtemp,t_half);
	t_full = t + tau;
	xtemp = x + tau*F3;
	F4 = derivsRK(xtemp,t_full);
	xout = x + tau/6.*(F1 + F4 + 2.*(F2+F3));
	return xout;
def rka(x,t,tau,err,derivsRK):

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
	tSave = t;  xSave = x;    # Save initial values
	safe1 = .9;  safe2 = 4.;  # Safety factors
	eps = np.spacing(1) # smallest value
	#%* Loop over maximum number of attempts to satisfy error bound
	maxTry = 100;
	for iTry in range (1,maxTry):
		#%* Take the two small time steps
		half_tau = 0.5 * tau;
		xTemp = rk4(xSave,tSave,half_tau,derivsRK);
		t = tSave + half_tau;
		xSmall = rk4(xTemp,t,half_tau,derivsRK);
		#%* Take the single big time step
		t = tSave + tau;
		xBig = rk4(xSave,tSave,tau,derivsRK);
		#%* Compute the estimated truncation error
		scale = err * (np.abs(xSmall) + np.abs(xBig))/2.;
		xDiff = xSmall - xBig;
		errorRatio = np.max( [np.abs(xDiff)/(scale + eps)] );
		#print safe1,tau,errorRatio
		#%* Estimate news tau value (including safety factors)
		tau_old = tau;
		tau = safe1*tau_old*errorRatio**(-0.20);
		tau = np.max([tau,tau_old/safe2]);
		tau = np.min([tau,safe2*tau_old]);
		#%* If error is acceptable, return computed values
		if errorRatio < 1 :
			xSmall = xSmall; # +  (xDiff)/15;
			return xSmall, t, tau
	#%* Issue error message if error bound never satisfied
	print('ERROR: Adaptive Runge-Kutta routine failed');
	return()
def gravrknbodies(s,t):
	#  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
	#  Inputs
	#    s      State vector [re(1) re(2) rj(1) rj(2) ve(1) ve(2) vj(1) vj(2)]
	#    t      Time (not used)
	#    mass(n>1)   Mass of objects in solar masses
	#  Output
	#    deriv  Derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
	# for an arbitrary number of bodies
	#* Compute acceleration
	# first get the number of bodies
	nbodies = len(mass); # one less since we have GM in the first point
	# now check to see if the length of s is 4*nbodies
	if len(s)/4. != nbodies:
		print(' number of bodies = #d ',nbodies)
		print(' length(s)/4 of bodies = ',len(s)/4.)
		print(' Wrong number of of bodies')
		return
	deriv = np.zeros(4*nbodies);
	# now cycle thru the bodies to get the acceleration
	for i in range (0,nbodies):
		for j in range(0,nbodies):
			if i !=j and  mass[i] > 0:  # only move if mass > 0
				xindex = 2*i;
				yindex = 2*i+1;
				xjndex = 2*j;
				yjndex = 2*j+1;
				vxindex = 2*(nbodies+i);
				vyindex = 2*(nbodies+i)+1;
				# must go 1 index more at stop
				deriv[vxindex:vyindex+1] = deriv[vxindex:vyindex+1] - \
				GM*mass[j]*(s[xindex:yindex+1]- s[xjndex:yjndex+1])/np.linalg.norm(s[xindex:yindex+1]-s[xjndex:yjndex+1])**3
				deriv[xindex:yindex+1] = s[vxindex:vyindex+1];
	#* Return derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
	return deriv

def rotate(th):
	# returns a rotation matrix for a given input angle "th"

	rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
	return rot

def lambdanwt(lam0, masses):
	# use Newton's method to find lambda for the Euler collinear solution
	itermax = 1000
	m1 = masses[0]
	m2 = masses[1]
	m3 = masses[2]
	lam = lam0 # set initial guess
	for it in range(0, itermax): # Newton's method
		lam_new = lam - (((m1+m2)*lam**5 + (3*m1+2*m2)*lam**4 + (3*m1+m2)*lam**3 - (m2+3*m3)*lam**2 - (2*m2+3*m3)*lam - (m2+m3)) /
			(5*(m1+m2)*lam**4 + 4*(3*m1+2*m2)*lam**3 + 3*(3*m1+m2)*lam**2 - 2*(m2+3*m3)*lam - (2*m2+3*m3)))
		lam = lam_new
	
	return lam

GM = 4*np.pi**2 # Grav. const. * Mass of Sun (au^3/yr^2)
option=int(input('Choose a solution:\n 1 Figure 8\n 2 Perturbed Figure 8\n 3 4 body braid\n \
	4 5 bodies\n 5 circle-o-planets\n 6 Lagrange\n 7 Euler Collinear: '))
if(option==1):
	# figure 8 - 3 bodies
	n_bodies=3
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	mass=np.zeros(n_bodies)
	r0[0,:] =[3.3030197, -0.82771837];  v0[0,:] = [1.587433767, 
	1.47221479];  mass[0] = 1;
	r0[1,:]= [-3.3030197, 0.82771837];  v0[1,:] = [1.587433767, 
	1.47221479];  mass[1] = 1;
	r0[2,:]= [0.0, 0.0];                v0[2,:] = [-3.174867535,-
	2.94442961]; mass[2] = 1;
elif(option==2):
	# figure 8 - 3 bodies - small perturbation
	n_bodies=3
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	mass=np.zeros(n_bodies)
	r0[0,:] = [3.3030197, -0.82771837]; v0[0,:] = [1.587433767, 
	1.47221479];  mass[0] = 1;
	r0[1,:]= [-3.3030197, 0.82771837];  v0[1,:] = [1.587433767, 
	1.47221479];  mass[1] = 1;
	r0[2,:]= [0.0, 0.1];                v0[2,:] = [-3.174867535,-
	2.94442961]; mass[2] = 1;
elif(option==3):
	# 4 bodies
	n_bodies=4
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	mass=np.zeros(n_bodies)
	r0[0,:] = [8.688745801755996, 0.0];  v0[0,:] = [0.0, 
	1.466058273458379 ];  mass[0] = 1;
	r0[1,:] = [-8.688745801755996, 0.0]; v0[1,:] = [0.0, -
	1.466058273458379];  mass[1] = 1;
	r0[2,:] = [0.0, 0.986648100464895];  v0[2,:] = [4.692245815597808, 
	0.0];   mass[2] = 1;
	r0[3,:] = [0.0, -0.986648100464895]; v0[3,:] = [-4.692245815597808, 
	0.0];  mass[3] = 1;
elif(option==4):
	# 5 bodies
	n_bodies=5
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	mass=np.zeros(n_bodies)
	r0[0,:] = [10.415422734260494, 0.0];                 v0[0,:] = [0.0, 1.488398966802320 ];                  mass[0] = 1;
	r0[1,:] = [2.763188135302640, 1.066362495707276];    v0[1,:] = [4.569045641873520 , -0.321470077002120];   mass[1] = 1;
	r0[2,:] = [2.763188135302640, -1.066362495707276];   v0[2,:] = [-4.569045641873520,  -0.321470077002120];  mass[2] = 1;
	r0[3,:] = [-7.970899502432867,   1.681699405991083]; v0[3,:] =[3.187338116419745,  -0.422729406399095];    mass[3] = 1;
	r0[4,:] = [-7.970899502432867,  -1.681699405991083]; v0[4,:] =[- 3.187338116419745,  -0.422729406399095];  mass[4] = 1;
elif(option==5):
	# circular orbit
	n_bodies = int(input('input number of bodies: '))
	if (n_bodies < 2):
		print('n_bodies must be greater than 1')
	r0 = np.zeros((n_bodies,2));
	v0 = np.zeros((n_bodies,2));
	mass = np.ones(n_bodies);
	# n objects in a circle
	vel0 = 0.;
	dtheta = 2*np.pi/(n_bodies);
	ang = np.linspace(0,2*np.pi-dtheta,n_bodies)
	# compute the speed in order to have a circular orbit
	for n in range(1,n_bodies):
		vel0 = vel0 + np.sin(ang[n]/2)/(1.-np.cos(ang[n]))/2;
	vel0 = np.sqrt(vel0*GM);
	period = 2*np.pi/vel0;
	print('Period = ',period,' years')
	for n in range(0,n_bodies):
		r0[n,:] = np.array([np.cos(ang[n]), np.sin(ang[n])]);
		v0[n,:] = vel0*np.array([-np.sin(ang[n]),np.cos(ang[n])]);
		mass[n] = 1;
elif(option==6):
	# 3 bodies in a triangle - Lagrange solution
	n_bodies=3
	mass = np.zeros(n_bodies)
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	speeds = np.zeros(n_bodies)
	accels = np.zeros(n_bodies)
	theta = (2*np.pi)/3

	# decide on initial condition
	s3 = np.array([1., 0.])
	mass[0] = 1
	mass[1] = 2
	mass[2] = 3
	mass_tot = np.sum(mass)

	# set initial positions and velocities
	r0[0,:] = ((mass[2]*rotate(-theta) - mass[1]*np.identity(2)).dot(s3)) / mass_tot
	r0[1,:] = ((mass[0]*np.identity(2) - mass[2]*rotate(theta)).dot(s3)) / mass_tot
	r0[2,:] = ((mass[1]*rotate(theta) - mass[0]*rotate(-theta)).dot(s3)) / mass_tot

	accels = np.arctan2(r0[:,1], r0[:,0])

	speeds[0] = 0.6*np.sqrt(GM)*(np.power((mass[1]**2 + mass[1]*mass[2]+ mass[2]**2), (3/4)) / mass_tot) * (1/np.sqrt(np.linalg.norm(r0[0,:])))
	speeds[1] = 0.6*np.sqrt(GM)*(np.power((mass[0]**2 + mass[0]*mass[2]+ mass[2]**2), (3/4)) / mass_tot) * (1/np.sqrt(np.linalg.norm(r0[1,:])))
	speeds[2] = 0.6*np.sqrt(GM)*(np.power((mass[0]**2 + mass[0]*mass[1]+ mass[1]**2), (3/4)) / mass_tot) * (1/np.sqrt(np.linalg.norm(r0[2,:])))

	v0[0,0] = -speeds[0]*np.sin(accels[0])
	v0[0,1] = speeds[0]*np.cos(accels[0])
	v0[1,0] = -speeds[1]*np.sin(accels[1])
	v0[1,1] = speeds[1]*np.cos(accels[1])
	v0[2,0] = -speeds[2]*np.sin(accels[2])
	v0[2,1] = speeds[2]*np.cos(accels[2])

else:
	# Euler collinear solution to  body problem
	n_bodies = 3
	mass = np.zeros(n_bodies)
	r0=np.zeros((n_bodies,2))
	v0=np.zeros((n_bodies,2))
	speeds = np.zeros(n_bodies)
	

	# decide on initial condition
	s3 = np.array([1., 0.])
	mass[0] = 1
	mass[1] = 2
	mass[2] = 3
	mass_tot = np.sum(mass)
	lam_guess = -2

	lambd = lambdanwt(lam_guess, mass)
	s1 = lambd*s3
	s2 = -(1+lambd)*s3

	# set initial position
	r0[0,:] = (mass[2]*s2 - mass[1]*s3)/mass_tot
	r0[1,:] = (mass[0]*s3 - mass[2]*s1)/mass_tot
	r0[2,:] = (mass[1]*s1 - mass[0]*s2)/mass_tot

	# set initial speeds
	speeds[0] = 0.6*np.sqrt(((mass[1]*GM/np.linalg.norm(r0[0,:]-r0[1,:])**2) +
	 (mass[2]*GM/np.linalg.norm(r0[0,:]-r0[2,:])**2))*np.linalg.norm(r0[0,:]))

	speeds[1] = 0.6*np.sqrt(((-mass[0]*GM/np.linalg.norm(r0[1,:]-r0[0,:])**2) +
	 (mass[2]*GM/np.linalg.norm(r0[1,:]-r0[2,:])**2))*np.linalg.norm(r0[1,:]))

	speeds[2] = 0.6*np.sqrt(((mass[0]*GM/np.linalg.norm(r0[2,:]-r0[0,:])**2) +
	 (mass[1]*GM/np.linalg.norm(r0[2,:]-r0[1,:])**2))*np.linalg.norm(r0[2,:]))

	accels = np.arctan2(r0[:,1], r0[:,0])

	# set initial velocities
	v0[0,0] = -speeds[0]*np.sin(accels[0])
	v0[0,1] = speeds[0]*np.cos(accels[0])
	v0[1,0] = -speeds[1]*np.sin(accels[1])
	v0[1,1] = speeds[1]*np.cos(accels[1])
	v0[2,0] = -speeds[2]*np.sin(accels[2])
	v0[2,1] = speeds[2]*np.cos(accels[2])



s0 = np.zeros((n_bodies,2))
vs0 = np.zeros((n_bodies,2))
for i in range(0,n_bodies):
	s0[i,:] = r0[i,:]*mass[i];
	vs0[i,:] = v0[i,:]*mass[i];
masst = np.sum(mass);
com=np.zeros((2,1))
vcom=np.zeros((2,1))
com = np.sum(s0[0:n_bodies,:],axis=0)/masst
vcom = np.sum(vs0[0:n_bodies,:],axis=0)/masst
r_com = np.linalg.norm(com); v_com=np.linalg.norm(vcom)
print('|Center of mass|=',r_com,' v_com=',v_com)
if(r_com > 1.0e-4 or v_com > 1.0e-4):
	print(' Warning, center of mass location or speed not zero')

adaptErr = 1.e-6; # Error parameter used by adaptive Runge-Kutta
time = 0.0
r = np.array([]);v=np.array([])
for n in range(0,n_bodies):
	r = np.append(r,r0[n,:])
	v = np.append(v,v0[n,:])
state = np.append(r,v);   # Used by R-K routines
#* Loop over desired number of steps using specified numerical method.
nStep = int(input('Enter number of steps: '))
tau = float(input('Enter time step (yr): '))
tmax = float(input('Enter the stop time (yr): '))
animation =int(input('Enter a 1 for animation: '));
method=1;
if method ==0:
	print(' Using rk4 method ');
elif method==1:
	print(' Using rka method, adaptErr = ',adaptErr);
# initialize arrays for plotting
xplot = np.array([]);   yplot=np.array([]); vxplot = np.array([]);   
vyplot=np.array([])
tplot=np.array([])
potential=np.array([]); kinetic=np.array([])
for iStep in range (0,nStep):
	#* Record position and energy for plotting.
	if(iStep ==0):
		xplot = np.append(xplot,r[0:2*n_bodies+1:2]);
		yplot = np.append(yplot,r[1:2*n_bodies+1:2]);
		vxplot = np.append(xplot,v[0:2*n_bodies+1:2]);
		vyplot = np.append(yplot,v[1:2*n_bodies+1:2]);
	else:
		xplot = np.vstack((xplot,r[0:2*n_bodies+1:2]));
		yplot = np.vstack((yplot,r[1:2*n_bodies+1:2]));
		vxplot = np.vstack((xplot,v[0:2*n_bodies+1:2]));
		vyplot = np.vstack((yplot,v[1:2*n_bodies+1:2]));
	tplot   = np.append(tplot,time)
	# potential is relative to the Sun's position
	kinetic_temp = 0.0;potential_temp = 0.0
	for n in range(0,n_bodies):
		for m in range(0,n_bodies):
			nxindex = 2*n;
			nyindex = 2*n+1;
			mxindex = 2*m;
			myindex = 2*m+1;
			if(nxindex != mxindex and n!= m):
				potential_temp =  potential_temp -  0.5 * GM * mass[n]*mass[m]/np.linalg.norm([r[nxindex]- r[mxindex], r[nyindex]-r[myindex]])
		# kinetic energy
		vxindex = 2*n;
		vyindex = 2*n+1;
		kinetic_temp = kinetic_temp+0.5*mass[n]*(v[vxindex]**2+v[vyindex]**2)
	potential=np.append(potential,potential_temp)
	kinetic=np.append(kinetic,kinetic_temp)
	if(time > tmax):
		print(' time exceeded tmax, breaking out of loop')
		break
	#* Calculate new position and velocity using desired method.
	# Note: In python, the matlab-type stop:end should be stop:end+1
	if method==0:
		state = rk4(state,time,tau,gravrknbodies);
		r[0:2*n_bodies] = state[0:2*n_bodies];   # 4th order Runge-Kutta
		v[0:2*n_bodies] = state[2*n_bodies:4*n_bodies];
		time = time + tau;
	elif method==1:
		[state,time,tau] = rka(state,time,tau,adaptErr,gravrknbodies);
		r[0:2*n_bodies] = state[0:2*n_bodies];   # adaptive 4th order Runge-Kutta
		v[0:2*n_bodies] = state[2*n_bodies:4*n_bodies];
# plotting routines
xmax = np.max(xplot);
xmin = np.min(xplot);
ymax = np.max(yplot);
ymin = np.min(yplot);

if (animation==1):
	# now animate
	plt.figure(4)
	plt.ion()
	plt.clf()
	for i in range(0,len(xplot)):
		plt.plot(xplot[0,:],yplot[0,:],'k+')
		plt.axis([xmin, xmax, ymin, ymax])
		plt.axis('equal');
		for n in range(0,n_bodies):
			plt.plot(xplot[1:i,n],yplot[1:i,n],'-');
			plt.plot(xplot[i,n],yplot[i,n],'ro');
			#        legend('m=3','m=4','m=5')
		plt.title('t = '+str(tplot[i])+' of '+str(tplot[-1]))
		plt.draw()
		plt.pause(0.01)
		if(i==len(xplot)-1):
			temp=input('Hit return to end:')
		else:
			plt.clf()
			# plt.ioff()
	# save the output figure if desired
	# plt.savefig('c4pA_ellip.png', bbox_inches='tight')
