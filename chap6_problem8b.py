import numpy as np
import matplotlib.pyplot as plt

# non Neumann demo - Richardson scheme
# Benjamin Klimko, PHYS 416, Spring 2018
# Program to plot the amplification factors as a function of phase angle kh

angle = np.arange(0,2*np.pi,np.pi/50.)
d=np.array([0.125, 0.25, 0.5, 0.6 ]) # d = Tau * kappa/h^2
r1=np.ones(len(angle))
# amplification factor as a function of d and angle
a1=np.zeros([len(angle),len(d)])
aa1=np.copy(a1)
x1=np.zeros([len(angle),len(d)])
y1=np.copy(x1)
a2=np.zeros([len(angle),len(d)])
aa2=np.copy(a2)
x2=np.zeros([len(angle),len(d)])
y2=np.copy(x2)
for j in range(0,len(d)):
	for i in range(0,len(angle)):
		coeffs = [1, 8*d[j]*np.sin(angle[i]/2)**2, -1]
		sols = np.roots(coeffs)
		a1[i,j]=sols[0]
		aa1[i,j]=np.abs(a1[i,j])
		a2[i,j] = sols[1]
		aa2[i,j] = np.abs(a2[i,j])
		#  for plotting
		x1[i,j] = a1[i,j]*np.cos(angle[i])
		y1[i,j] = a1[i,j]*np.sin(angle[i])
		x2[i,j] = a2[i,j]*np.cos(angle[i])
		y2[i,j] = a2[i,j]*np.sin(angle[i])
print(str(d[0:]))
plt.figure(1)
plt.clf()
plt.plot(x1[:,0],y1[:,0],x1[:,1],y1[:,1],x1[:,2],y1[:,2],x1[:,3],y1[:,3], x2[:,0],y2[:,0],x2[:,1],y2[:,1],x2[:,2],y2[:,2],x2[:,3],y2[:,3])
leg=['d='+str(d[0]),' d='+str(d[1]),' d='+str(d[2]),' d='+str(d[3])]
plt.grid(True)
plt.legend(leg)
plt.axis('equal')
plt.plot(np.cos(angle),np.sin(angle),'r--',linewidth=2)
plt.grid(True)
plt.title('A')
plt.figure(2)
plt.clf()
plt.plot(angle,a1[:,0],angle,a1[:,1],angle,a1[:,2],angle,a1[:,3], angle,a2[:,0],angle,a2[:,1],angle,a2[:,2],angle,a2[:,3])
plt.grid(True)
plt.ylabel(' A')
plt.xlabel('angle kh (radians)')
plt.legend(leg)
plt.show()