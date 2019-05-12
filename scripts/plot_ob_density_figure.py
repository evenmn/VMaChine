import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

numberOfDimensions = 2

files = ["../data/int1/onebody/VMC/2D/12P/0.100000w/SGD_MC1048576.dat",
         "../data/int1/onebody/VMC/2D/12P/0.280000w/SGD_MC1048576.dat",
         "../data/int1/onebody/VMC/2D/12P/0.500000w/SGD_MC1048576.dat",
         "../data/int1/onebody/VMC/2D/12P/1.000000w/SGD_MC1048576.dat",
         ]
         
label = ["$\omega=0.1$", 
         "$\omega=0.28$", 
         "$\omega=0.5$", 
         "$\omega=1.0$",
         ]
         
maxRadius = [20,10,10,10]

sns.set()

for i in range(len(files)):
    data = np.loadtxt(files[i])
    data = data[-100*maxRadius[i]:]
    r = np.linspace(0,maxRadius[i],len(data))
    data /= (r**(numberOfDimensions-1))
    data /= np.sum(np.nan_to_num(data))
    data2 = data[::-1]
    data2 = np.concatenate((data2,data))
    r2 = np.linspace(-maxRadius[i],maxRadius[i],2*len(data))
    plt.plot(r2, data2, '.', markersize=3, alpha=1.0, label=label[i])

plt.xticks([])
plt.yticks([])
plt.grid()
plt.show()
