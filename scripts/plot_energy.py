import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from scipy.signal import savgol_filter

asymptote = 0.073839

files = ["../data/int1/quantumdot/energy/VMC/2D/2P/0.010000w/ADAM_MC1048576.dat", 
         "../data/int1/quantumdot/energy/RBM/2D/2P/0.010000w/ADAM_MC1048576.dat",
         "../data/int1/quantumdot/energy/RBMSJ/2D/2P/0.010000w/ADAM_MC1048576.dat",
         "../data/int1/quantumdot/energy/RBMPJ/2D/2P/0.010000w/ADAM_MC1048576.dat"
         ]

label = ["VMC", 
         "RBM", 
         "RBM+SJ",
         "RBM+PJ"
         ]

for i in range(len(files)):
    data = np.loadtxt(files[i])
    x = np.linspace(0, len(data) - 1, len(data))
    x = savgol_filter(x, 51, 2)
    plt.plot(x, data, label=label[i])

label_size = {"size":"14"}

plt.axhline(asymptote, linestyle='--', color='r', label="Exact")

plt.xlabel("Iteration",**label_size)
plt.ylabel("Dimensionless energy",**label_size)
plt.legend()
plt.grid()
plt.show()
'''
fig, ax = plt.subplots()
ax.plot(x, data1, label="Brute-force Metropolis")
ax.plot(x, data2, label="Metropolis-Hastings")
ax.plot(x, data3, label="Gibbs' sampling")

ax.axhline(asymptote, linestyle='--', color='r', label="Exact")
plt.xlabel("Iteration")
plt.ylabel("Energy")

plt.xlabel(u'Iteration', fontname = 'Times New Roman', size = 14)
plt.ylabel(u'Energy [a.u.]', fontname = 'Times New Roman', size = 14)


#legend = plt.legend(loc='upper right', shadow=False, fontsize='large')
axins = zoomed_inset_axes(ax, 10, loc=10)

axins.axhline(asymptote, linestyle='--', color='r', label="Exact")
axins.plot(x, data1)
axins.plot(x, data2)
axins.plot(x, data3)

x1, x2, y1, y2 = 960, 1001, 2.96, 3.1

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
#axins.grid()
ax.grid()
ax.legend(loc='upper right')

plt.show()
'''
