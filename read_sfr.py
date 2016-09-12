import matplotlib.pyplot as plt
import numpy as np

sfr = np.loadtxt('sfr_history.dat')
time = np.zeros(sfr.shape[0],dtype = int)

for i in range(0,sfr.shape[0]):
	time[i] = i

n_stars = np.round(sfr/0.35)
plt.plot(time/1000, sfr)
plt.ylabel("SFR (M$_\odot$ Myr$^{-1}$)")
plt.xlabel("Time (Gyr)")
plt.xlim([0,time[sfr.shape[0]-1]/1000])
plt.show()

plt.clf()

sn_rate = np.loadtxt('sn_rate.dat')

plt.plot(time, sn_rate,marker='.', ls='None')
plt.ylabel("SN Rate (Myr$^{-1}$)")
plt.xlabel("Time (Gyr)")
plt.xlim([0,time[sfr.shape[0]-1]])
plt.show()