import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


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
    r_mat = np.zeros((numBins, numBins))
    for j in range(numBins):
        for k in range(numBins):
            r_mat[j,k] = (j*j + k*k)**((numberOfDimensions-1)/2)

    data /= r_mat
    data /= np.sum(np.nan_to_num(data))
    return data
    
    
def plot(data, radius):
    plt.figure()
    #sns.heatmap(data, cmap="YlGnBu")
    plt.imshow(data, cmap=plt.cm.jet, extent=[0,radius,radius,0])
    plt.colorbar()
    plt.xlabel("Radial distance from particle 1")
    plt.ylabel("Radial distance from particle 2")


def main():
    maxRadius = [10]
    newRadius = [4]
    
    methods   = ['VMC']
    dims      = ['2']
    particles = ['2']
    omegas    = ['0.280000']
    
    for method in methods:
        for dim in dims:
            for particle in particles:
                i = 0
                for omega in omegas:
                    fileName = generateFileName(method, dim, particle, omega)
                    data = np.loadtxt(fileName)
                    data = crop(data, maxRadius[i], newRadius[i])
                    data = norm(data, int(dim))
                    plot(data, newRadius[i])
                    saveFigure(method, dim, particle, omega)
                    i += 1
    plt.show()
        

if __name__ == '__main__':
    main()
