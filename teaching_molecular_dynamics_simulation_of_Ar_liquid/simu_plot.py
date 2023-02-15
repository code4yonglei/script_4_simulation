#!/usr/bin/env python3

import matplotlib.pylab as plt
import numpy as np

def moving_avg(x):
    y = np.zeros(len(x))
    current_sum=0.0
    for i in range(len(x)):
        current_sum+= x[i]
        y[i] = current_sum/(i+1.0)
    return y

data=np.loadtxt('simu_ener.txt')
t = data[:,1]
T = data[:,2]
P = data[:,3]
PE= data[:,4]
KE= data[:,5]
TE= data[:,6]
px= data[:,7]
py= data[:,8]
pz= data[:,9]


# Energy
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t, PE, '-b', lw=2.0, label='PE')
plt.plot(t, KE, '-r', lw=2.0, label='KE')
plt.plot(t, TE, '-g', lw=2.0, label='TE')
plt.title('Argon Liquid Simulation (Energy)', fontsize=14)
plt.ylabel('Energy', fontsize=12)
# Energy Drift
plt.subplot(2,1,2)
TE_cavg = moving_avg(data[10:,6])
drift = (data[10:,6] - TE_cavg)/256.0
plt.plot(data[10:,1], drift, '-k', lw=2.0, label='Energy Drift')
drift_avg = moving_avg(drift)
print(drift_avg)
plt.plot(data[10:,1], drift_avg, '--r', lw=2.0, label='Avg Energy Drift')
plt.xlabel('Time [units]', fontsize=12)
plt.ylabel('Energy Drift per atom', fontsize=12)
plt.legend(loc='best', fontsize=12)
plt.savefig('simu_Energy.png',dpi=300)
plt.show()


# Temperature
plt.figure(2)
t_sample=data[10:,1]
T_sample=data[10:,2]
T_cavg=moving_avg(T_sample)
plt.plot(t, T, '-b', lw=2.0, label='Instant Temperature')
plt.plot(t_sample, T_cavg, '-r', lw=2.0, label='Averaged Temperature')
plt.legend(loc='best', fontsize=12)
plt.title('Argon Liquid Simulation (Temperature)', fontsize=14)
plt.xlabel('Time [units]', fontsize=12)
plt.ylabel('Temperature [K]', fontsize=12)
plt.savefig('simu_Temperature.png', dpi=300)
plt.show()


# Pressure
plt.figure(3)
t_sample=data[10:,1]
P_sample=data[10:,3]
P_cavg=moving_avg(P_sample)
print(P_cavg)
plt.plot(t, P, '-b', lw=2.0, label='Instant Pressure')
plt.plot(t_sample, P_cavg, '-r', lw=2.0, label='Averaged Pressure')
plt.legend(loc='best', fontsize=12)
plt.title('Argon Liquid Simulation (Pressure)', fontsize=14)
plt.xlabel('Time [units]',fontsize=12)
plt.ylabel('Pressure [MPa]',fontsize=12)
plt.savefig('simu_Pressure.png',dpi=300)
plt.show()


# Momentum
plt.figure(4)
plt.plot(t, px, '-b', lw=2.0, label='Px')
plt.plot(t, py, 'or', lw=2.0, label='Py')
plt.plot(t, pz, '-g', lw=2.0, label='Pz')
plt.legend(loc='best', fontsize=12)
plt.title('Argon Liquid Simulation (Momentum)', fontsize=14)
plt.xlabel('Time [units]', fontsize=12)
plt.ylabel('Momentum', fontsize=12)
plt.savefig('simu_Momentum.png',dpi=300)
plt.show()
