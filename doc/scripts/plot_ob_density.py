import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def exact(r1, w):
    '''Exact solution without interaction for arbitrary w'''
    return (2*r1+1)*np.exp(- w * r1**2)

# Numerical values
data1 = np.loadtxt("../data/onebody_VMC_P6_D2_w1.000000_MC16777216.dat")
data2 = np.loadtxt("../data/onebody_VMC_P6_D2_w1.000000_MC16777216.dat")

label_size = {"size":"16"}

data1 = data1[-500:]
data2 = data2[-500:]
r1 = np.linspace(0, 5, len(data1))

exact1 = exact(r1, w=1.0)/np.sum(exact(r1, w=1.0))

data1/=(r1)
data2/=(r1)

data1  = data1/np.sum(np.nan_to_num(data1))
data2  = data2/np.sum(np.nan_to_num(data2))

plt.plot(r1, exact1, '--r', linewidth=1.0, label="Exact")
plt.plot(r1, data1, '-.', markersize=2, alpha=1.0, label="VMC")
plt.plot(r1, data2, ':', markersize=2, alpha=1.0, label="NQS")
#plt.title("Two non-interacting particles with $\omega=0.5$")
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best", fontsize=16)
plt.grid()
#plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
