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
    fileName  = "../data/int0/"
    fileName += system   + "/"
    fileName += "twobody" + "/"
    fileName += method   + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += "GD_MC65536.dat"
    return fileName


def saveFigure(system, method, dim, particle, omega):
    fileName  = "../plots/int0/"
    #fileName += system   + "/"
    fileName += "twobody" + "/"
    #fileName += "VMC" + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += method   + "_GD_MC2pow30.png"
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
    
def norm(data, numberOfDimensions, numberOfParticles):
    numBins = len(data)
    v = radial(numBins, numberOfDimensions)
    xx, yy = np.meshgrid(v, v, sparse=True)

    data /= np.multiply(xx,yy)
    data /= np.sum(np.nan_to_num(data))
    
    data = np.where(data > 0.000016, 0, data)
    
    data *= 4e4
    
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
    
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def plot(data, radius):
    size = 24
    size_ticks = 16
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams.update({'figure.autolayout': True})

    fig, ax = plt.subplots(figsize=(8,6))
                 
    img = ax.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    cbar = fig.colorbar(img, fraction=0.046, pad=0.04)#, format=ticker.FuncFormatter(fmt))
    cbar.set_label(r'$\rho(r_i,r_j)$', rotation=270, labelpad=40, y=0.45, **label_size)
    cbar.ax.tick_params(labelsize=size_ticks)
    
    plt.tight_layout()
    
    ax.set_xlabel("$r_j$", **label_size)
    ax.set_ylabel("$r_i$", **label_size)
    ax.tick_params(labelsize=size_ticks)
    
    ticks = [-3,-2,-1, 0, 1, 2, 3]
        
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    plt.grid()



def main():
    maxRadius = [3]
    newRadius = [3]

    systems   = ['quantumdot']
    methods   = ['VMC']
    dims      = ['2']
    particles = ['2']
    omegas    = ['1.000000']      

    i = 0
    for system in systems:
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
                        plot(data, newRadius[0])
                        #saveFigure(system, method, dim, particle, omega)
                        i += 1
    plt.show()


if __name__ == '__main__':
    main()

