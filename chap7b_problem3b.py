import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# Benjamin Klimko, PHYS 416, Spring 2018
# traffic simulation program using particles
# speed of the cars depends on distance to the car in front
# discretized car simulation, cars move to the right and hit a wall

xmax = 400      # length of system
car_len = 10    # length of each car
vmax = 25.0     # speed limit
tmax = 20       # time limit
dt = 0.2        # time step
itmax = int(tmax/dt) # number of iterations
# define a small car shape
carx = np.array([-car_len/2.,  car_len/2.,  car_len/2.,car_len/3.,  car_len/4., -car_len/4., -  car_len/3., -car_len/2., -car_len/2.])

cary = tmax/100.*np.array([0., 0., car_len/3., car_len/3., 2*car_len/3., 2*car_len/3., car_len/3.,  car_len/3., 0.])
# intialize car positions
# xcar = 1:2*car_len:xmax
xcar = np.arange(-xmax/4.,0,car_len)
x0 = np.copy(xcar)
number_of_cars = len(xcar) # number of cars
# set the initial density
density=np.zeros(number_of_cars)
speed=np.zeros(number_of_cars)
density[0:-1] = car_len/(xcar[1:number_of_cars]-xcar[0:number_of_cars-1])
density[-1] = car_len/(xmax/2.-xcar[-2])/2.
# set the initial speed
speed[0:] = np.maximum(0.,vmax*(1.-density[0:]))
xcar_time = xcar
density_time=density
speed_time=speed
time=np.array([0.0])
# plotting
plt.ion()
plt.figure(1)
plt.plot(xcar[0]+carx,cary,'b')
plt.axis([0.,xmax, 0., tmax+car_len])
plt.xlabel('x')
plt.ylabel('time')
for n in range(0,number_of_cars):
	plt.plot(xcar[n]+carx,cary,'b')
# now loop through time
xcar_time = xcar
xcar_analytic = np.copy(xcar_time)
density_time = density
speed_time = speed
for i in range(1,itmax):
	xcar = xcar + dt*speed
	time = np.append(time,i*dt)
	# density
	density[0:-1] = car_len/(xcar[1:number_of_cars]-xcar[0:number_of_cars-
	1])
	density[-1] = car_len/(xmax/2-xcar[number_of_cars-1])/2.
	density=np.minimum(1.,density)
	# speed - keep it positive
	speed[0:] = np.maximum(0,vmax*(1-density[0:]))
	xcar_time   =np.vstack((xcar_time,xcar))
	temp = np.zeros(number_of_cars)
	for car in range(0,number_of_cars):
		if (i*dt) < (-x0[car]/vmax):
			temp[car] = x0[car]
		else:
			temp[car] = vmax*(i*dt) - 2*np.sqrt(-x0[car]*vmax*i*dt)

	xcar_analytic = np.vstack((xcar_analytic, temp))
	density_time=np.vstack((density_time,density))
	speed_time  =np.vstack((speed_time,speed))
	plt.clf()
	plt.plot(xcar[0]+carx,cary+i*dt,'b')
	plt.axis([-xmax/2., xmax/2., 0., tmax+car_len])
	plt.xlabel('x')
	plt.ylabel('time')
	for n in range(1,number_of_cars):
		plt.plot(xcar[n]+carx,cary+i*dt,'b')
	plt.draw()
	# plotting the car locations
	plt.plot(xcar_time[:,0],time,'r')
	plt.plot(xcar_analytic[:,0], time, 'g')
	for n in range(1,number_of_cars):
		plt.plot(xcar_time[:,n],time,'r')
		plt.plot(xcar_analytic[:,n], time, 'g')
	plt.draw()
	plt.pause(dt)
	# some IDEs need this to work (eg, PyCharm)
	# if i == itmax-1:
	# 	temp=input('Hit any key to stop: ')

# plt.savefig('c7bq3b.png', bbox_inches='tight') # save figure if desired