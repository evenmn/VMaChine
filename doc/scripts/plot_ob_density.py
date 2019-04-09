import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def exact(r1, w):
    '''Exact solution without interaction for arbitrary w'''
    return (2*r1+1)*np.exp(- w * r1**2)


files = ["../data/onebody_VMC_SGD_P8_D3_INT1_w1.000000_MC4194304.dat",
         "../data/onebody_NQS_SGD_P8_D3_INT1_w1.000000_MC4194304.dat",
         "../data/onebody_VMC_SGD_P8_D3_INT1_w0.500000_MC4194304.dat",
         "../data/onebody_NQS_SGD_P8_D3_INT1_w0.500000_MC4194304.dat",
         "../data/onebody_VMC_SGD_P8_D3_INT1_w0.280000_MC4194304.dat",
         "../data/onebody_NQS_SGD_P8_D3_INT1_w0.280000_MC4194304.dat",
         "../data/onebody_VMC_SGD_P8_D3_INT1_w0.100000_MC4194304.dat",
         "../data/onebody_NQS_SGD_P8_D3_INT1_w0.100000_MC4194304.dat"
         ]
         
label = ["VMC, $\omega=1.0$", 
         "RBM, $\omega=1.0$", 
         "VMC, $\omega=0.5$", 
         "RBM, $\omega=0.5$", 
         "VMC, $\omega=0.28$", 
         "RBM, $\omega=0.28$", 
         "VMC, $\omega=0.1$", 
         "RBM, $\omega=0.1$"
         ]

for i in range(len(files)):
    data = np.loadtxt(files[i])
    data = data[-500:]
    r = np.linspace(0,5,len(data))
    data /= (r*r)
    data /= np.sum(np.nan_to_num(data))
    plt.plot(r, data, '.', markersize=2, alpha=1.0, label=label[i])


label_size = {"size":"16"}

#exact1 = exact(r1, w=1.0)/np.sum(exact(r1, w=1.0))
#plt.plot(r1, exact1, '--r', linewidth=1.0, label="Exact")

#plt.title("Two non-interacting particles with $\omega=0.5$")
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
plt.legend(loc="best", fontsize=16, markerscale=5)
plt.grid()
#plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
