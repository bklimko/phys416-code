import numpy as np
import matplotlib.pyplot as plt
"""
interp - Program to interpolate data using Lagrange
polynomial to fit quadratic to three data points
clear all; help interp;   Clear memory and print header
Initialize the data points to be fit by quadratic
"""
def intrpf(xi,x,y):
	"""
	Function to interpolate between data points
	using Lagrange polynomial (quadratic)
	Inputs
	x    Vector of x coordinates of data points (3 values)
	y    Vector of y coordinates of data points (3 values)
	xi   The x value where interpolation is computed
	Output
	yi   The interpolation polynomial evaluated at xi
	"""
	# Calculate yi = p(xi) using Lagrange polynomial
	yi = (xi-x[1])*(xi-x[2])/((x[0]-x[1])*(x[0]-x[2]))*y[0]\
	+ (xi-x[0])*(xi-x[2])/((x[1]-x[0])*(x[1]-x[2]))*y[1] \
	+ (xi-x[0])*(xi-x[1])/((x[2]-x[0])*(x[2]-x[1]))*y[2];
	return (yi);
x = [0.,0.,0.];
y = [0.,0.,0.]
print ('Enter data points as x,y pairs (e.g., 1,2)');
for i in range(0,3):
	temp1,temp2 =  input("Enter data points: ").split(',');
	x[i] = float(temp1);
	y[i] = float(temp2);
# Establish the range of interpolation (from x_min, x_max)
temp1,temp2 = input("Enter range of x values as x_min, x_max: ").split(',');
xr=np.array([float(temp1),float(temp2)]);
# Find yi for the desired interpolation values xi using
# the function intrpf
nplot = 100;     #% Number of points for interpolation curve
xi= np.linspace(xr[0], xr[1], nplot,endpoint=True)
yi= np.linspace(0, 0, nplot,endpoint=True)
for i in range(0,nplot):
	xi[i] = xr[0] + (xr[1]-xr[0])*(i-1)/(nplot-1);
	yi[i] = intrpf(xi[i],x,y);  #% Use intrpf function to interpolate
# Plot the curve given by (xi,yi) and mark original data points
plt.figure(1);plt.clf() # open a window and clear it
plt.plot(x,y,'*',label="Data points")
plt.plot(xi,yi,'-',label="Interpolation");
plt.xlabel('x');
plt.ylabel('y');
plt.title('Three point interpolation');
plt.legend()
plt.show()