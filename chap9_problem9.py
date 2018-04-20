import numpy as np
import matplotlib.pyplot as plt 
# Benjamin Klimko, PHYS 416, Spring 2018
# Program to find the spectral radius of a matrix involved with the Lax scheme for the advection equation
def power1(A,X,eps,max1):
	# function to compute the maximum eigenvalue (lambda) of a Matrix A and
	# its corresponding eigen vector (eigenvector)
	# X is some base vector row matrix, input
	# eps is the tolerance you want the eigenvalue to be computed to
	# max1 is the maximum number of iteractions allowed
	lamda=0.0
	cnt=1
	xnew=X
	xold=0*X
	err=np.linalg.norm(X)
	xnew = xnew/max(abs(xnew))
	while err>eps and cnt < max1:
		xold=xnew
		ck = max(abs(xnew))
		xnew=np.dot(A,xnew)/ck
		cnt = cnt+1
		err = np.linalg.norm(xold-xnew)
	if (cnt >=max1):
		print('max number of iterations exceeded')
	eigenvector = xnew/np.linalg.norm(xnew)
	lamda = ck
	return lamda,eigenvector

# set Lax scheme matrices
B = np.zeros((51,51))
C = np.zeros((51,51))
x = np.random.randint(5, size=(51,1))
th = np.linspace(0.1,1.5, 15)
powval = np.array([])
eigenvals = np.array([])
onenorm = np.array([])
infnorm = np.array([])

for idx in range(1,50):
	B[idx, idx-1] = -1
	B[idx, idx+1] = 1
	C[idx, idx-1] = 1
	C[idx, idx+1] = 1

# matrix edge cases
B[0,1] = 1
B[0, -1] = -1
B[-1,0] = 1
B[-1,-2] = -1
C[0,1] = 1
C[0, -1] = 1
C[-1,0] = 1
C[-1,-2] = 1

for frac in th:
	print('tau/h =', str(frac))
	A = ((1/2)*C) - ((frac/2)*B)
	# x=np.array([[3.], [2.], [1.]])
	
	eigenvalues = np.linalg.eigvals(A)

	# get max eigenvalue
	# emax = np.amax(eigenvalues)
	emax = eigenvalues.max()
	eigenvals = np.append(eigenvals, emax)

	eigenvalue1,eigenvector1=power1(A,x,1.0e-3,20) # power method
	powval = np.append(powval, eigenvalue1)

	eigenvalue2 = np.linalg.norm(A, 1) # 1-norm
	onenorm = np.append(onenorm, eigenvalue2)

	eigenvalue3 = np.linalg.norm(A, np.inf) # infinity norm
	infnorm = np.append(infnorm, eigenvalue3)

	print('max eigenvalue from numpy eig =',emax)

	print('max eigenvalue from powermethod =',eigenvalue1)

	print('max eigenvalue from 1-norm =',eigenvalue2)

	print('max eigenvalue from inf-norm =',eigenvalue3)

# plot results
plt.figure()
plt.plot(th, powval, label='Power method')
plt.plot(th, eigenvals, label='Numpy eig function')
plt.plot(th, onenorm, label='1-norm')
plt.plot(th, infnorm, label='Infinity norm')
plt.legend()
plt.show()