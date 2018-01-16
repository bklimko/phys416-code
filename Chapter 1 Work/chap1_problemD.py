from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Original program by Dr. Frank Toffoletto, edited by Benjamin Klimko, PHYS 416 Spring 2018
The program has a user definable input "levels", which defines the plotted shape by how much
recursion happens in the move0 function
This program also calculates the area and perimeter of the von Koch curve and plots them for various levels of recursion

"""

def move0(x,y,angle,s, level, pm):
	# routine to add a side of length s at angle angle to arrays x and y
	if level == 0:
		# calculate as normal
		x = np.append(x, x[-1]+s*np.cos(np.pi/180.*angle))
		y = np.append(y, y[-1]+s*np.sin(np.pi/180.*angle))
		pm += s

	else:
		# recursively define the mapping F -> F - F ++ F - F as defined in the problem handout
		[x,y, pm] = move0(x,y,angle,s/3, level-1, pm)
		angle = angle - rotation
		[x,y, pm] = move0(x,y, angle, s/3, level-1, pm)
		angle = angle + rotation
		angle = angle + rotation
		[x,y, pm] = move0(x,y,angle,s/3, level-1, pm)
		angle = angle - rotation
		[x,y, pm] = move0(x,y,angle, s/3, level-1, pm)
	return x,y, pm

# arrays to hold perimeter and area info for plotting purposes
p = []
a = []
# origin point for triangle method of finding area
x1 = 0
y1 = 0

# loop through various levels of recursion in order to find how area and perimeter vary with level
for level in range(0,8):   # level = 0 gives triangle, = 1 gives star, >1 gives snowflake
	x=0; y=0; angle=0; #% starting point and direction
	n_sides = 3
	rotation = 180/n_sides; #% rotation angle
	s = 1; #% unit of step length
	# starts off at zero
	x=[0.];
	y=[0.];
	perimeter = 0  # initialize perimeter and area for each run
	area  = 0

	for i in range (0,n_sides):
		
		[x,y, perimeter] =move0(x,y,angle,s, level, perimeter); #% move one unit
		angle = angle + rotation; #% rotation in degrees
		angle = angle + rotation; #% rotation in degrees

	p = np.append(p, perimeter)
	
	# calculate area for each level
	for tri in range(0, len(x)-2):
		mat = np.array([[x1, y1, 1], [x[tri], y[tri], 1], [x[tri+1], y[tri+1], 1]])
		area += 0.5*np.linalg.det(mat)

	a = np.append(a, abs(area))
	
	
#% now plot the result

plt.subplot(2,1,1)
plt.plot(range(0,8), p)
plt.xlabel('Iterations')
plt.ylabel('Perimeter')
plt.title('Iterations vs Perimeter')

plt.subplot(2, 1, 2)
plt.plot(range(0,8), a)
plt.xlabel('Iterations')
plt.ylabel('Area')
plt.title('Iterations vs Area')
plt.tight_layout()
plt.show()