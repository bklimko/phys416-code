import numpy as np 
import matplotlib.pyplot as plt 

# Benjamin Klimko, PHYS 416, Spring 2018
# Program to compute equilibrium position for two dimensional spring-mass system

# set constants and initial conditions


g = 9.81 # gravity
m = 0.1 # mass 
D = 0.1
L1 = 0.1
L2 = 0.1
Dmat = np.zeros((2,2))

thres = 1e-5 # threshold for convergence of Newton's method

x0 = 4
y0 = 3
x = x0
y = y0
xvec = np.array([x,y])

imax  = 1000

for step in range(0, imax):
	# Create Jacobian
	Dmat[0,0] = ((20*(np.sqrt((x-0.1)**2 + y**2)-0.1))/np.sqrt((x-0.1)**2 + y**2) + 10*(np.sqrt(x**2 + y**2)-0.1)/np.sqrt(x**2 + y**2)
		+10*x**2/(x**2+y**2) + 5*(2*x -0.2)**2/((x-0.1)**2 + y**2) - 10*x**2*(np.sqrt(x**2 + y**2)-0.1)/((x**2+y**2)**(3/2))
		- 5*(np.sqrt((x-0.1)**2+y**2)-0.1)*(2*x-0.2)**2/(((x-0.1)**2+y**2)**(3/2)))

	Dmat[0,1] = ((10*x*y)/(x**2 + y**2) + (10*y*(2*x- 0.2))/((x-0.1)**2 + y**2) - (10*y*(np.sqrt((x-0.1)**2 + y**2)-0.1)*(2*x-0.2))/ 
		(((x-0.1)**2 + y**2)**(3/2)) - (10*x*y*(np.sqrt((x-0.1)**2 + y**2)-0.1))/((x**2 + y**2)**(3/2)))

	Dmat[1,0] = ((10*x*y)/(x**2 + y**2) + (10*y*(2*x- 0.2))/((x-0.1)**2 + y**2) - (10*y*(np.sqrt((x-0.1)**2 + y**2)-0.1)*(2*x-0.2))/ 
		(((x-0.1)**2 + y**2)**(3/2)) - (10*x*y*(np.sqrt((x-0.1)**2 + y**2)-0.1))/((x**2 + y**2)**(3/2)))

	Dmat[1,1] = (20*y**2/((x-0.1)**2+y**2) + 20*(np.sqrt((x-0.1)**2+y**2)-0.1)/np.sqrt((x-0.1)**2 + y**2)+ 
		10*(np.sqrt(x**2 + y**2)-0.1)/np.sqrt(x**2 + y**2) + 10*y**2/(x**2 + y**2) - 20*y**2*(np.sqrt((x-0.1)**2+y**2)-0.1)/((x-0.1)**2+y**2)**(3/2)
		- 10*y**2*(np.sqrt(x**2 + y**2)-0.1)/(x**2+y**2)**(3/2))

	# compute forces
	Fx = 10*x*(np.sqrt(x**2 + y**2)-0.1)/np.sqrt(x**2 + y**2) + 10*(2*x -0.2)*(np.sqrt((x-0.1)**2 + y**2)-0.1)/np.sqrt((x-0.1)**2 + y**2)
	Fy = 10*y*(np.sqrt(x**2 + y**2)-0.1)/np.sqrt(x**2 + y**2) + 20*y*(np.sqrt((x-0.1)**2 + y**2)-0.1)/np.sqrt((x-0.1)**2 + y**2) - m*g 
	gradient = np.array([Fx, Fy])
	# Newton's method
	xvec = xvec - gradient.dot(np.linalg.inv(Dmat))
	x = xvec[0]
	y = xvec[1]

# graph equilibrium position of system
plt.figure()
plt.plot(x,y, 'ko')
plt.plot([0, x], [2*y, y])
plt.plot([.1, x], [2*y, y])
plt.plot([-.5, .5],[2*y, 2*y], 'k')
plt.title('Equilibrium position: x = '+ str(np.around(x, decimals=5))+ ', y = '+ str(np.around(y, decimals=5)))
plt.show()