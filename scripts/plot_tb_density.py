import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def radial(size, numberOfDimensions):
    n = np.arange(size)
    if numberOfDimensions == 2:
        return 2*n+1
    elif numberOfDimensions == 3:
        return 3*n*(n+1)+1
    else:
        return np.nan

def generateFileName(method, dim, particle, omega):
    fileName  = "../data/int1/twobody/"
    fileName += method   + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += "ADAM_MC1048576.dat"
    return fileName


def saveFigure(method, dim, particle, omega):
    fileName  = "../plots/int1/twobody/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += method   + "_ADAM_MC2pow28.png"
    plt.savefig(fileName)


def crop(data, maxRadius, newRadius):
    length = len(data)
    newLength = int(length * (newRadius / maxRadius))
    newData = data[0:newLength, 0:newLength]
    return newData

    
    
def norm(data, numberOfDimensions):
    numBins = len(data)
    v = radial(numBins, numberOfDimensions)
    xx, yy = np.meshgrid(v, v, sparse=True)
    data /= xx
    data /= yy
    data /= np.sum(np.nan_to_num(data))
    return data
    
    
def rotate(data):
    data1 = np.rot90(data,0)
    data2 = np.rot90(data,-1)
    data3 = np.rot90(data,-2)
    data4 = np.rot90(data,-3)
    
    dataTop = np.c_[data2,data1]
    dataBottom = np.c_[data3,data4]
    data = np.concatenate((dataBottom,dataTop))
    return data
    


def plot(data, radius):
    size = 16
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.gcf().subplots_adjust(left=0.18)

    plt.figure()
    #sns.heatmap(data, cmap="YlGnBu")                                           
    plt.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    plt.colorbar()
    plt.xlabel("Radial distance from particle $j$", **label_size)
    plt.ylabel("Radial distance from particle $i$", **label_size)



def main():
    maxRadius = [40]
    newRadius = [6]

    methods   = ['RBMPJ']
    dims      = ['2']
    particles = ['42']
    omegas    = ['1.000000'] #'1.000000','0.500000','0.280000','0.100000']      

    i=0
    for method in methods:
        for dim in dims:
            for particle in particles:
                for omega in omegas:
                    fileName = generateFileName(method, dim, particle, omega)
                    data = np.loadtxt(fileName)
                    data = crop(data, maxRadius[i], newRadius[i])
                    data = norm(data, int(dim))
                    data = rotate(data)
                    plot(data, newRadius[i])
                    saveFigure(method, dim, particle, omega)
                    i += 1
    plt.show()


if __name__ == '__main__':
    main()

