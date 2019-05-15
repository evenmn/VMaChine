import numpy as np
import matplotlib.pyplot as plt

time_vmc = [4.7485, 8.551862, 16.69836, 34.40414, 68.92126, 125.4514, 214.0844, 353.5306]
time_rbm = [3.4,6.551862, 13.69836, 30.40414, 65.92126, 119.4514, 207.0844, 345.5306]
orbitals = np.arange(1,9)

label_size = {'size':'14'}
plt.plot(orbitals, time_rbm, label='RBM')
plt.plot(orbitals, time_vmc, '--', label='VMC')
plt.xlabel('Number of orbitals',**label_size)
plt.ylabel('CPU-time [s]',**label_size)
plt.legend(loc='best', fontsize=14)
plt.grid()
plt.show()
