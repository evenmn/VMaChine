import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def exact(r1):
    return np.exp(-r1**2)

# Numerical values
data1 = np.loadtxt("../data/OB.dat")
#data2 = np.loadtxt("../data/OB_w_interaction.dat")

data1 = data1[9500:]

label_size = {"size":"14"}

r1 = np.linspace(0, 5, len(data1))

data2 = exact(r1)
data1 = data1/np.sum(data1)
data2 = 1.65*data2/np.sum(data2)

data1/=(r1)
plt.plot(r1, data1, '.', markersize=2, label="Interaction")
plt.plot(r1, data2, '--', linewidth=1.0, label="W")
plt.title("2 particles with interaction")
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
#plt.legend(loc="best")
plt.grid()
plt.savefig("../html/ob_2^22.png")
plt.show()
