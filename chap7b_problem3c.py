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

# now loop through time
xcar_time = xcar
xcar_analytic = np.copy(xcar_time)
density_time = density
speed_time = speed
times = np.zeros(number_of_cars)
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
	
	for car in range(0,number_of_cars):
		if xcar_time[-1,car] > 0 and xcar_time[-2,car] < 0:
			times[car] = (i*dt)- dt/2
		

	density_time=np.vstack((density_time,density))
	speed_time  =np.vstack((speed_time,speed))


plt.figure()
plt.plot(x0, times, label='Numerical solution')
xtime = (-4*x0)/vmax
plt.plot(x0, xtime, label='Analytic solution')
plt.legend()
plt.title('Starting position vs Time to reach intersection')
plt.xlabel('Starting position (arb length)')
plt.ylabel('Time (s)')
plt.show()
temp=input('Hit any key to stop: ')
