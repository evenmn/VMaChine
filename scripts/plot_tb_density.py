import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import simps

plt.style.use("bmh")

def radial(size, numberOfDimensions):
    n = np.arange(size)
    if numberOfDimensions == 2:
        return 2*n+1
    elif numberOfDimensions == 3:
        return 3*n*(n+1)+1
    else:
        return np.nan

def generateFileName(system, method, dim, particle, omega):
    fileName  = "../data/int1/"
    fileName += system   + "/"
    fileName += "twobody" + "/"
    fileName += method   + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += "ADAM_MC1048576.dat"
    return fileName


def saveFigure(system, method, dim, particle, omega):
    fileName  = "../plots/int1/"
    #fileName += system   + "/"
    fileName += "twobody" + "/"
    #fileName += "VMC" + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += method   + "_ADAM_MC1048576.png"
    print(fileName)
    plt.savefig(fileName)

def crop(data, maxRadius, newRadius):
    length = len(data)
    newLength = int(length * (newRadius / maxRadius))
    newData = data[0:newLength, 0:newLength]
    return newData

def remove_cross(data):
    lengthHalf = int(len(data)/2)
    data[lengthHalf] = data[lengthHalf-1]
    data[lengthHalf+1] = data[lengthHalf+2]
    data[:,lengthHalf] = data[:,lengthHalf-1]
    data[:,lengthHalf+1] = data[:,lengthHalf+2]
    return data
    
def norm(data, numberOfDimensions, numberOfParticles, factor=1e5):
    numBins = len(data)
    v = radial(numBins, numberOfDimensions)
    xx, yy = np.meshgrid(v, v, sparse=True)
    data /= np.multiply(xx,yy)
    #data /= np.sum(np.nan_to_num(data))
    norm_const = simps(simps(data, v), v)
    data /= norm_const
    data *= factor
    return data * numberOfParticles
    
def rotate(data):
    data = data[2:,2:]

    data1 = np.rot90(data,0)
    data2 = np.rot90(data,-1)
    data3 = np.rot90(data,-2)
    data4 = np.rot90(data,-3)
    
    dataTop = np.c_[data2,data1]
    dataBottom = np.c_[data3,data4]
    data = np.concatenate((dataBottom,dataTop))
    return data
    
def cut(data, threshold=0.00001):
    return np.where(data > threshold, 0, data)
    
def ticks(radius):
    if radius % 2 == 0 or radius < 2:
        return [-radius, -radius/2., 0, radius/2., radius]
    elif radius % 3 == 0:
        return [-radius, -2*(radius/3.), -radius/3., 0, radius/3., 2*(radius/3.), radius]
    elif radius == 5:
        return [-5, -3, -1, 1, 3, 5]
    elif radius == 25:
        return [-25, -15, -5, 5, 15, 25]
    elif radius == 35:
        return [-35, -17.5, 0, 17.5, 35]
    else:
        print("Warning: Ticks out of bounds")
        return []
    
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def plot(data, radius):
    size = 28
    size_ticks = 20
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams.update({'figure.autolayout': True})

    fig, ax = plt.subplots(figsize=(8,6))
                 
    img = ax.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    cbar = fig.colorbar(img, fraction=0.046, pad=0.04)#, format=ticker.FuncFormatter(fmt))
    cbar.set_label(r'$\rho(r_i,r_j)$', rotation=90, labelpad=10, y=0.5, **label_size)
    cbar.ax.tick_params(labelsize=size_ticks)
    
    plt.tight_layout()
    
    ax.set_xlabel("$r_j$", **label_size)
    ax.set_ylabel("$r_i$", **label_size)
    ax.tick_params(labelsize=size_ticks)
    
    tick = ticks(radius)
        
    ax.set_xticks(tick)
    ax.set_yticks(tick)
    plt.grid()
    #plt.show()


if __name__ == '__main__':
    system   = 'quantumdot'

    maxRadius = [55]
    newRadius = [50]

    methods   = ['RBMSJ']#,'RBM']#,'RBMSJ','RBMPJ']
    dims      = ['2']
    particles = ['6']
    omegas    = ['0.010000']      

    i = 0
    for method in methods:
        for dim in dims:
            for particle in particles:
                for omega in omegas:
                    fileName = generateFileName(system, method, dim, particle, omega)
                    data = np.loadtxt(fileName)
                    data = crop(data, maxRadius[0], newRadius[0])
                    data = norm(data, int(dim), int(particle))
                    data = rotate(data)
                    data = remove_cross(data)
                    #data = cut(data, 0.06)
                    plot(data, newRadius[0])
                    saveFigure(system, method, dim, particle, omega)
                    plt.show()
                    i += 1

