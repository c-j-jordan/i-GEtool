import numpy as np
import matplotlib.pyplot as plt
import math


fname = 'star_masses1.dat'
data = None
file = open(fname, 'r')
#fileread = file.read()
#file.close()
#stripped = fileread.strip(" ")
#data = stripped.split()
#stripped = None
#fileread = None
#np.asarray(data)



data = file.read()

splitted = data.split()

masses = np.asarray(splitted)

mass = np.zeros(np.size(masses),dtype=float)

logmass = np.zeros(np.size(masses),dtype=float)

# Reads in z y x

for i in range(0,np.size(masses)):
            mass[i] = masses[i]

file.close()
