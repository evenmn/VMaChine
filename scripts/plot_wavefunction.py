import numpy as np
import matplotlib.pyplot as plt

def generateFileName(system, method, dim, particle, omega):
    fileName  = "../data/int1/"
    fileName += system   + "/"
    fileName += "weights" + "/"
    fileName += method   + "/"
    fileName += dim      + "D/"
    fileName += particle + "P/"
    fileName += omega    + "w/"
    fileName += "ADAM_MC1048576.dat"
    return fileName
    
def evaluate_wavefunction(data, method, particle, dim, omega, positions):
    if method == VMC:
        gauss = np.exp(-0.5 * omega * data[0,0] * np.inner(positions, positions))
        jastrow = 

def main():
    maxRadius = [15]
    newRadius = [15]

    systems   = ['quantumdot']
    methods   = ['RBM']
    dims      = ['2']
    particles = ['6']
    omegas    = ['0.100000']      

    i=0
    for system in systems:
        for method in methods:
            for dim in dims:
                for particle in particles:
                    for omega in omegas:
                        fileName = generateFileName(system, method, dim, particle, omega)
                        data = np.loadtxt(fileName)

                        plot(data, newRadius[0])
                        saveFigure(system, method, dim, particle, omega)
                        i += 1
    plt.show()


if __name__ == '__main__':
    main()
