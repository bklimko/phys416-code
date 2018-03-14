import numpy as np 
import matplotlib.pyplot as plt 

#Benjamin Klimko, PHYS 416, Spring 2018
# Program that plots total length of a mass-spring system acting as a model of a short polymer molecule

def equipos(ratio):
	#output the matrix and vector for the equilib. position equation
	mat = np.array([[-2*ratio - 1, ratio, 1],[ratio, -2*ratio - 1, ratio],[1, ratio, -ratio - 1]])
	vec = np.array([[1],[-1], [-ratio -1]])
	return mat, vec 


ratio = np.arange(1e-3, 1e3) # define ratio k1/k2 over a range of values
length = np.zeros(len(ratio))

# for each k1/k2 ratio record the total length (i.e. the position x3 since the first mass is fixed at x0 = 0)
for q in range(0, len(ratio)):
	[m, v] = equipos(ratio[q])
	# matrix inverse times vector gives position vector
	length[q] = (np.linalg.inv(m).dot(v))[2]

# plot results
plt.figure()
plt.semilogx(ratio, length)
plt.xlabel(r'$k_1$/$k_2$')
plt.ylabel('Total System length')
plt.show()
