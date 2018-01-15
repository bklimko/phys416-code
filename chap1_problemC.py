from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Original program by Dr. Frank Toffoletto, edited by Benjamin Klimko, PHYS 416 Spring 2018
The program has a user definable input "levels", which defines the plotted shape by how much
recursion happens in the move0 function 

"""

def move0(x,y,angle,s, level):
	# routine to add a side of length s at angle angle to arrays x and y
	if level == 0:
		# calculate as normal
		x = np.append(x, x[-1]+s*np.cos(np.pi/180.*angle))
		y = np.append(y, y[-1]+s*np.sin(np.pi/180.*angle))

	else:
		# recursively define the mapping F -> F - F ++ F - F as defined in the problem handout
		[x,y] = move0(x,y,angle,s/3, level-1)
		angle = angle - rotation
		[x,y] = move0(x,y, angle, s/3, level-1)
		angle = angle + rotation
		angle = angle + rotation
		[x,y] = move0(x,y,angle,s/3, level-1)
		angle = angle - rotation
		[x,y] = move0(x,y,angle, s/3, level-1)
	return x,y


x=0; y=0; angle=0; #% starting point and direction
level = 3 #% level = 0 gives triangle, = 1 gives star, >1 gives snowflake
n_sides = 3
rotation = 180/n_sides; #% rotation angle
s = 1; #% unit of step length
# starts off at zero
x=[0.];
y=[0.];

for i in range (0,n_sides):
	
	[x,y]=move0(x,y,angle,s, level); #% move one unit
	angle = angle + rotation; #% rotation in degrees
	angle = angle + rotation; #% rotation in degrees
	
#% now plot the result

plt.figure(1),plt.clf()
plt.plot(x,y)
plt.axis([-1.5*s, 1.5*s,-1.5*s,1.5*s])
plt.axis('equal')
#py.axis('tight')
plt.axis('off')
plt.show()