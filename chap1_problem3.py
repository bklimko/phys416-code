import numpy as np
import matplotlib.pyplot as plt  # import statements for plotting purposes

"""
Benjamin Klimko, PHYS 416 Spring 2018
This program recreates a series of plots from the textbook as an exercise to become comfortable with Python
Running the program will produce four subplots on the same figure; the non-labeled functions are as follows:
--dashed red line: exp(-x/4)
--polar plot: abs(sin(3x))
"""

# create the linspaces of the coordinates for the various functions
xcoords1 = np.linspace(0, 19, 1000)  # for graphed functions 1, 3, and 4
xcoords2 = np.linspace(1, 19, 10)  # for graphed function 2
theta = np.linspace(0, 2*np.pi, 500)


# apply the functions to the coordinates from above
func1 = np.exp(-xcoords1/4)*np.sin(xcoords1)
func2 = np.exp(-xcoords2/4)*np.sin(xcoords2)
func3 = np.exp(-xcoords1/4)
func4 = np.exp(-xcoords1)*(np.sin(xcoords1))**2
func5 = np.abs(np.sin(3*theta))

# graph the first subplot
plt.subplot(2, 2, 1)
plt.plot(xcoords1, func1)
# set axis labels and title
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('f(x)=exp(-x/4)*sin(x)')
# set axis tick labels and grid
plt.xticks([0, 5, 10, 15, 20], ('0', '5', '10', '15', '20'))
plt.yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1], ('-0.4', '-0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'))
plt.grid(True)

# graph the second subplot
plt.subplot(2, 2, 2)
plt.plot(xcoords1, func1)
plt.plot(xcoords2, func2, 'go', markerfacecolor='none')  # markerfacecolor lets you make the green circles hollow
plt.plot(xcoords1, func3, 'r--')
# give axis tick labels, axis labels, title, and grid
plt.xticks([0, 5, 10, 15, 20], ('0', '5', '10', '15', '20'))
plt.yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1], ('-0.4', '-0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'))
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('f(x)=exp(-x/4)*sin(x)')
plt.grid(True)

# graph the third subplot
plt.subplot(2, 2, 3)
plt.semilogy(xcoords1, func4, '.')
# set axis tick labels, axis labels, title, and grid
plt.xticks([0, 5, 10, 15], ('0', '5', '10', '15'))
plt.yticks([10e-10, 10e-5, 1], (r'$10^{-10}$', r'$10^{-5}$', r'$10^0$'))
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(r'f(x)=exp(-x)$*sin(x)^2$')
plt.grid(True)
# TODO FIGURE OUT HOW TO MAKE THIS LOOK MORE RIGHT

# graph the fourth subplot, give axis labels, set axis ticks
plt.subplot(2, 2, 4, projection='polar')
plt.plot(theta, func5)
ax = plt.gca()

# TODO stuck with editing the details of the polar plot

# display the graphs
plt.show()
