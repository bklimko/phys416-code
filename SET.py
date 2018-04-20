import numpy as np 
import matplotlib.pyplot as plt 

#SET- program to model a single electron transistor and plot its I-V characteristics
# Benjamin Klimko, PHYS 416, Spring 2018

def tunnel_rate(delF, temp):
	# function to determine the tunneling rate for a given change in system free energy and temperature
	rate = delF/(e**2 * R * (1-np.exp(-delF/(kb*temp))))
	return rate 

# define parameters-- can change to allow for user inputs

Cd = 0.06e-18 # Drain capacitance 
Cs = np.copy(Cd) # Source capacitance
Cg = 0.23e-18 # Gate capacitance
Cb = np.copy(Cg) # Backgate capacitance

Vg = 0 # Gate voltage
Vb = 0 # Backgate voltage
R = 1e6 # Junction resistance 
T = 300 # Absolute temperature
n = 3 # number of states
kb = 1.381e-23 # Boltzmann constant (J/K)

# uncomment below to allow for user input of parameters instead of hardcoding

# Cd = float(input('Enter Cd in Farads: ')) # Drain capacitance 
# Cs = float(input('Enter Cs in Farads: ')) # Source capacitance
# Cg = float(input('Enter Cg in Farads: ')) # Gate capacitance
# Cb = float(input('Enter Cb in Farads: ')) # Backgate capacitance
# Vg = float(input('Enter Vg in Volts: ')) # Gate voltage
# Vb = float(input('Enter Vb in Volts: ')) # Backgate voltage
# R = float(input('Enter R in Ohms: ')) # Junction resistance 
# T = float(input('Enter T in degrees Kelvin: ')) # Absolute temperature

e = 1.6e-19 # elementary charge
C = Cd + Cs + Cg + Cb # total system capacitance 
Vdd = np.linspace(-1, 1) # different drain voltages for testing
Ids = np.zeros(len(Vdd)) # matrix to store drain-source current for plotting
idx = 0 # indexing variable for loop

for volt in Vdd: # MAIN LOOP-- go through every drain voltage
	delta_F = (e/C) *(e/2 + volt*(C-Cd) - Vg*Cg - Vb*Cb - n*e) # calculate change in free energy related to tunneling event
	prob = ((1/C)*(Vg*Cg + Vb*Cb + volt*Cd + n*e))/3 # occupation probability as in Ismail and Abdelrassoul

	Ids[idx] = 2*e*prob*tunnel_rate(delta_F, T) # calculate Ids as in Ismail and Abdelrassoul
	idx += 1 # update index variable


# plot results
plt.plot(Vdd, Ids)
plt.xlabel(r'$V_{dd} (V)$')
plt.ylabel(r'$I_{ds} (A)$')
plt.title('I-V Characteristics of SET')
plt.show()

