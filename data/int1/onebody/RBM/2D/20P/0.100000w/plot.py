import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("../0.280000w/SGD_MC1048576.dat")
data = data[-1000:]
r = np.linspace(0,10,len(data))
data /= r
data /= np.sum(np.nan_to_num(data))
plt.plot(r, data, '.', markersize=1, alpha=1.0, label=r"$\omega=0.28$")
plt.show()

