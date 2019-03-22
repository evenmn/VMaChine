import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

asymptote = 20.19
'''
data = np.loadtxt("../data/energy.txt")
x = np.linspace(0, len(data) - 1, len(data))

plt.plot(x, data, label="Calculated")
plt.axhline(asymptote, linestyle='--', color='r', label="Exact")
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.grid()
plt.legend()
#plt.axis([-10, 210, 1.5, 6])
#plt.savefig('figure.png')
plt.show()

'''
#data1 = np.loadtxt("../data/energy_NQS_SGD_P2_D2_INT1_w0.280000_MC4194304.dat")
#data2 = np.loadtxt("../data/energy_NQS_SGD_P2_D2_INT1_w1.000000_MC33554432.dat")
#data3 = np.loadtxt("../data/energy_VMC_SGD_P2_D2_INT1_w0.500000_MC262144.dat")
data4 = np.loadtxt("../data/energy_VMC_SGD_P6_D2_INT1_w1.000000_MC262144.dat")

x = np.linspace(0, len(data4) - 1, len(data4))

label_size = {"size":"14"}

plt.axhline(asymptote, linestyle='--', color='r', label="Exact")
#plt.plot(x, data1, label="VMC")
#plt.plot(x, data2, label="NQS")
#plt.plot(x, data3, label="VMC")
plt.plot(x, data4, label="NQS")

plt.xlabel("Iteration",**label_size)
plt.ylabel("Energy [a.u.]",**label_size)
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
