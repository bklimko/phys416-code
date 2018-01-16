import numpy as np
import matplotlib.pyplot as plt 

"""
Benjamin Klimko, PHYS 416 Spring 2018
The program has a user definable input hardcoded (N) and produces plots of starting value vs computed values and starting value
vs iterations for the hailstone (3n+1) problem as performed in problem A of this problem set. 
"""

# initialize the start vector s, computed number vector f,  iteration vector g, input vector nums, and counter count
N = 200
nums = np.arange(1, N)
s = []
f = []
g = np.zeros(N-1)

# for each num in the input vector
for idx in nums:
	count = 0
	#code borrowed from my chap1_problemA.py 
	# make a copy of idx to keep track of what the starting number is
	n = int(np.copy(idx))


	# as long as n is not 1 continue the algorithm
	while n != 1:
		# add the current value of n to vector f and the value of idx to vector s for each time through the while loop
		s.append(idx)
		f.append(n)
		# check if n is even or odd-- if even divide by 2; if odd multiply by 3 and add 1
		if n % 2 == 0:
			n /= 2
		else:
			n = (3*n) + 1

		# increment the count for each time through the while loop
		count += 1

	s.append(idx)
	f.append(n)
	g[idx-1] = count

# plot starting value vs computed value and starting value vs iterations
plt.subplot(2, 1, 1)
plt.plot(s, f, '.')
strfn = lambda x: str(x)
vecfn = np.vectorize(strfn)
plt.xticks(np.arange(1,N+1, 10), vecfn(np.arange(1,N+1, 10)))
ax = plt.gca()
ax.set_xlim(1, N)
ax.set_ylim(0, max(f))
plt.tight_layout()

plt.subplot(2, 1, 2)
plt.plot(nums, g)
plt.xticks(np.arange(1,N,20), vecfn(np.arange(1,N,20)))
plt.yticks(np.arange(0,N,20), vecfn(np.arange(0,N,20)))
ax2 = plt.gca()
ax2.set_xlim(1, N)
ax2.set_ylim(0, max(g))

plt.show()