import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
#sns.set()

numberOfDimensions = 3

def radial(size):
    n = np.arange(size)
    if numberOfDimensions == 2:
        return 2*n+1
    elif numberOfDimensions == 3:
        return 3*n*(n+1)+1
    else:
        return np.nan


def exact(r1, w):
    '''Exact solution without interaction for given w'''
    return (2*r1*0+1)*np.exp(- w * r1**2)

files = ["../data/int0/onebody/VMC/3D/2P/1.000000w/GD_MC1048576.dat",
         "../data/int0/onebody/RBM/3D/2P/1.000000w/GD_MC1048576.dat",
         #"../data/int1/quantumdot/onebody/RBMSJ/3D/20P/1.000000w/ADAM_MC1048576.dat",
         #"../data/int1/quantumdot/onebody/RBMPJ/3D/20P/1.000000w/ADAM_MC1048576.dat",
         ]
         
label = ["VMC",
         "RBM",
         "RBM+SJ",
         "RBM+PJ"
         ]
         
line_style = ["-",
              "--", 
              "-.", 
              ":"
              ]
         
maxRadius = [10,
             10,
             10,
             10
             ]

limit = [0.0032, 
         0.0032, 
         0.0030, 
         0.0030
         ]

for i in range(len(files)):
    data = np.loadtxt(files[i])
    data /= radial(len(data))
    data /= np.sum(np.nan_to_num(data))
    data /= maxRadius[i]
    data /= np.sum(data)
    #data *= int(len(data)/1000)
    r = np.linspace(0,maxRadius[i],len(data))
    #data = np.where(data > 1.0005, 0, data) 
    data[:np.argmax(data)] = np.where(data[:np.argmax(data)] < limit[i], 0, data[:np.argmax(data)])
    indices = np.where(data == 0)[0]
    data = np.delete(data, indices)
    r = np.delete(r, indices)
    data /= np.sum(data)
    #r = r[np.argmax(data):]
    #data = data[np.argmax(data):]
    
    #data = savgol_filter(data, 101, 2)
    
    plt.plot(r, data, line_style[i], markersize=1, label=label[i])

size = 16
label_size = {"size":str(size)}
plt.rcParams["font.family"] = "Serif"
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.18)
#plt.tight_layout()

exact1 = exact(r, w=1.0)/np.sum(exact(r, w=1.0))
plt.plot(r, exact1, '--r', linewidth=1.0, label="Exact")

plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$(r)", **label_size)
plt.legend(loc="best", fontsize=size)
plt.axis([-0.1,5.1,-0.0001,0.0045])
plt.grid()
#plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
