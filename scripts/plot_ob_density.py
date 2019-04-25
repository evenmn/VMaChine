import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

numberOfDimensions = 2

def exact(r1, w):
    '''Exact solution without interaction for given w'''
    return (2*r1+1)*np.exp(- w * r1**2)

files = ["../data/int1/onebody/RBM/2D/20P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/20P/1.000000w/SGD_MC1048576.dat",
         "../data/int1/onebody/VMC/2D/20P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/6P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBM/2D/2P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/12P/0.500000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/6P/1.000000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/2P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBM/2D/20P/0.500000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/RBM/2D/20P/1.000000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/VMC/2D/20P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/VMC/2D/2P/0.100000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/VMC/2D/20P/0.500000w/SGD_MC1048576.dat",
         #"../data/int1/onebody/VMC/2D/20P/1.000000w/SGD_MC1048576.dat",
         ]
         
label = ["RBM",
         #"RBMPJ",
         "VMC",
         #"RBMPJ, $\omega=0.1$",
         #"RBMPJ, $\omega=0.28$", 
         #"RBMPJ, $\omega=0.5$", 
         #"RBMPJ, $\omega=1.0$", 
         #"RBM, $\omega=0.28$",
         #"VMC, $\omega=0.1$", 
         #"VMC, $\omega=0.28$", 
         #"VMC, $\omega=0.5$", 
         #"VMC, $\omega=1.0$",
         ]
         
maxRadius = [25,25,25]

for i in range(len(files)):
    data = np.loadtxt(files[i])
    data = data[-100*maxRadius[i]:]
    r = np.linspace(0,maxRadius[i],len(data))
    data /= (r**(numberOfDimensions-1))
    data /= np.sum(np.nan_to_num(data))
    plt.plot(r, data, '.', markersize=1, alpha=1.0, label=label[i])

label_size = {"size":"16"}
plt.rcParams["font.family"] = "serif"
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.18)
#plt.tight_layout()

#exact1 = exact(r1, w=1.0)/np.sum(exact(r1, w=1.0))
#plt.plot(r1, exact1, '--r', linewidth=1.0, label="Exact")

#plt.title("$\omega=0.1$", **label_size)
plt.xlabel("r / r$_0$", **label_size)
plt.ylabel(r"$\rho$(r) / r$_0^{\,-2}$", **label_size)
plt.legend(loc="best", fontsize=16, markerscale=5)
plt.grid()
#plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
