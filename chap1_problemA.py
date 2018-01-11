import numpy as np 
import matplotlib.pyplot as plt 

# set input value here and initialize iteration counter
n = 27
original_n = np.copy(n)
count = 0

# create arrays to hold both iteration values and computed numbers
iteration = np.array([0])
computed = np.array([n])

# as long as n is not 1 continue the algorithm
while n != 1:
	# check if n is even or odd-- if even divide by 2; if odd multiply by 3 and add 1
	if n%2 == 0:
		n /= 2
	else:
		n = (3*n) + 1

	# increment counter, add to iteration array, and add the updated value of n to computed array
	count += 1
	iteration = np.append(iteration, count)
	computed = np.append(computed, n)


# once n has become 1 plot the iterations vs computed number, add title and axis labels
plt.plot(iteration, computed)
plt.xlabel('Iterations')
plt.ylabel('Computed Number')
plt.title('Iterations vs Computed Number for N = {0}'.format(str(original_n)))

plt.show()

