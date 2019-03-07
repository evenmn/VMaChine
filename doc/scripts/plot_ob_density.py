import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def exact(r1):
    return np.exp(-r1**2)

# Numerical values
data1 = np.loadtxt("../data/onebody_PJ_P2_D2_MC1048576.dat")
data2 = np.loadtxt("../data/onebody_NQS_P2_D2_MC1048576.dat")

label_size = {"size":"14"}

data1 = data1[99501:]
data2 = data2[99501:]
r1 = np.linspace(0, 5, len(data1))

exact1 = exact(r1)/np.sum(exact(r1))
data1  = data1/np.sum(data1)
data2  = data2/np.sum(data2)

data1/=(r1)
plt.plot(r1, data1, '.', markersize=2, label="PJ")
plt.plot(r1, data2, '.', markersize=2, label="NQS")
plt.plot(r1, exact1, '--', linewidth=1.0, label="W")
plt.title("2 particles with interaction")
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best")
plt.grid()
plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
