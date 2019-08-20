import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set()

numberOfDimensions = 2

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
    return (2*r1+1)*np.exp(- w * r1**2)

files = ["../data/int1/onebody/VMC/2D/20P/0.280000w/ADAM_MC1048576.dat",
         #"../data/int1/onebody/RBM/2D/20P/0.280000w/ADAM_MC1048576.dat",
         "../data/int1/onebody/RBMSJ/2D/20P/0.280000w/ADAM_MC1048576.dat",
         #"../data/int1/onebody/RBMPJ/2D/20P/0.280000w/ADAM_MC1048576.dat",
         ]
         
label = ["VMC",
         #"RBM",
         "RBM+SJ",
         "RBM+PJ"
         ]
         
line_style = ["-",
              "--", 
              "-.", 
              ":"
              ]
         
maxRadius = [30,
             30,
             30,
             30
             ]

limit = [0.00013, 
         0.000125, 
         0.000125, 
         0.000
         ]

for i in range(len(files)):
    data = np.loadtxt(files[i])
    data /= radial(len(data))
    data /= np.sum(np.nan_to_num(data))
    data /= maxRadius[i]
    data *= int(len(data)/1000)
    r = np.linspace(0,maxRadius[i],len(data))
    #data /= np.sum(data)
    data = np.where(data > 1.0005, 0, data) 
    data[:np.argmax(data)] = np.where(data[:np.argmax(data)] < limit[i], 0, data[:np.argmax(data)])
    indices = np.where(data == 0)[0]
    data = np.delete(data, indices)
    r = np.delete(r, indices)
    #r = r[np.argmax(data):]
    #data = data[np.argmax(data):]
    
    plt.plot(r, data, line_style[i], markersize=1, label=label[i])

size = 16
label_size = {"size":str(size)}
plt.rcParams["font.family"] = "Serif"
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.18)
#plt.tight_layout()

#exact1 = exact(r1, w=1.0)/np.sum(exact(r1, w=1.0))
#plt.plot(r1, exact1, '--r', linewidth=1.0, label="Exact")

plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$(r)", **label_size)
plt.legend(loc="best", fontsize=size)
plt.grid()
#plt.savefig("../html/onebody_PJ_NQS_P2_D2_MC1048576.png")
plt.show()
