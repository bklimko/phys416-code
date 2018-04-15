import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
# relaxa - Program to solve the Laplace equation using
# the SOR method on a square grid
# this version makes an animation as it iterates
# It also reverses the coordinate system so that the first index is y and the
# second is x, which is in line with the grid structure
# Benjamin Klimko, PHYS 416, Spring 2018
#* Initialize parameters (system size, grid spacing, etc.)

N = int(input('Enter number of grid points on a side: '))
angle = np.linspace(0,2*np.pi,N)
hphi = (2*np.pi)/N 
L = 2          # System size (r goes to 2)
hr = L/(N-1)    # Grid spacing
r = np.arange(0,N)*hr   # r coordinate

xx,yy = np.meshgrid(r,angle) # for plotting
plot_interval =50 # interval to plot animation, setting it smaller slows the program down alot
animate=int(input(' Enter a 1 to animate while iterating:'))
#* Select over-relaxation factor 
omegaOpt = 2./(1.+np.sin(np.pi/N))  # Theoretical optimum
print('Theoretical optimum omega = ',omegaOpt)
omega = float(input('Enter desired omega: '))
#* Set initial guess as first term in separation of variables soln.

phi = np.zeros((N,N)) # row = angle, column = radius value
analyt_phi = np.zeros((N,N))

#* Set boundary conditions
# first index is the row and second index is column
phi[:,-1] = 0 # outer circle held at 0
idx = 0 
cosan = np.cos(angle)
for ang in cosan:
	phi[idx, int(N/2)] = ang 
	analyt_phi[idx, int(N/2)] = ang 
	idx += 1

mask = np.ones((N,N))
# set areas that do not get recalculated (hold outer box and inner box at respective potentials)
mask[:,0:int(N/2)] = 0
mask[:,-1] = 0
#* Loop until desired fractional change per iteration is obtained

iterMax = N**2          # Set max to avoid excessively long runs
changeDesired = 1.0e-4   # Stop when the change is given fraction
print('Desired fractional change = ',changeDesired)
change = np.array([])
for iter in range(0,iterMax):
	changeSum = 0.0
	## SOR method ##
	for i in range(1,N):        # Loop over interior points only
		for j in range(0,N):
			if mask[j,i] == 1: # check if the point should have its solution computed
				# if j != 99:
				# 	term1 = (phi[j,i-1]+phi[j,i+1])/(hr**2)
				# 	term2 = (1/(i*hr))*(phi[j,i+1]-phi[j,i-1])/(2*hr)
				# 	term3 = (1/(i*hr))**2 * (phi[j-1,i] - phi[j+1, i])/(hphi**2)
				# 	denom = (2/hr**2 + 2/((i*hr)**2 * hphi))
				# 	newphi = omega*((term1+term2+term3)/denom) +  (1-omega)*phi[j,i]
				# 	changeSum = changeSum + abs(1-phi[j,i]/newphi)
				# 	phi[j,i] = newphi
				# else:
				# 	term1 = (phi[j,i-1]+phi[j,i+1])/(hr**2)
				# 	term2 = (1/(i*hr))*(phi[j,i+1]-phi[j,i-1])/(2*hr)
				# 	term3 = (1/(i*hr))**2 * (phi[j-1,i] - phi[0, i])/(hphi**2)
				# 	denom = (2/hr**2 + 2/((i*hr)**2 * hphi))
				# 	newphi = omega*((term1+term2+term3)/denom) +  (1-omega)*phi[j,i]
				# 	changeSum = changeSum + abs(1-phi[j,i]/newphi)
				# 	phi[j,i] = newphi
				phi[j,i] = (4/(3*hr*i) - (i*hr)/3)*np.cos(hphi*j)
	#* Check if fractional change is small enough to halt the iteration
	# change = np.append(change,changeSum/(N-2)**2)
	if( iter > 1 ):
		# print('After %d iterations, fractional change = %f'%( iter,change[-1]))
		break 
	# if( change[-1] < changeDesired ):
	# 	print('Desired accuracy achieved after %d iterations'%iter)
	# 	print('Breaking out of main loop')
	# 	break
	# animate
	if(animate ==1 and iter%plot_interval<1):
		fig = plt.figure(2)   # Clear figure 2 window and bring forward
		plt.clf()
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(xx, yy, phi, rstride=1, cstride=1, 
		cmap=cm.jet,linewidth=0, antialiased=False)
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('potential after '+str(iter)+' iterations')
		plt.draw()
		plt.show()
		plt.pause(0.1)

#* Plot final estimate of potential as contour and surface plots

fig = plt.figure(2)   # Clear figure 2 window and bring forward
plt.clf()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, yy, phi, rstride=1, cstride=1, 
cmap=cm.jet,linewidth=0, antialiased=False)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('potential after '+str(iter)+' iterations')
plt.show()
