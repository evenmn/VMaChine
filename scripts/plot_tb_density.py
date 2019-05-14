import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# Load file
fileNames = ["../data/int1/twobody/VMC/2D/2P/0.100000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/2P/0.280000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/2P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/2P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/6P/0.100000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/6P/0.280000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/6P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/6P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/12P/0.100000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/12P/0.280000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/12P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/12P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/20P/0.100000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/20P/0.280000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/20P/0.500000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/20P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/30P/0.100000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/30P/0.280000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/30P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/30P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/42P/0.100000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/42P/0.280000w/ADAM_MC1048576.dat",
             #"../data/int1/twobody/VMC/2D/42P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/42P/1.000000w/ADAM_MC1048576.dat",
             
             "../data/int1/twobody/VMC/2D/56P/0.100000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/56P/0.280000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/56P/0.500000w/ADAM_MC1048576.dat",
             "../data/int1/twobody/VMC/2D/56P/1.000000w/ADAM_MC1048576.dat",
            ]


maxRadius = [10, 10, 10, 10,
             15, 15, 15, 15,
             25, 25,
             30,
             35, 35, 35, 35,
             40, 40, 40,
             45, 45, 45, 45,
             50, 50, 50, 50,
            ]


for i in range(len(fileNames)):
    print(i)
    data = np.loadtxt(fileNames[i])

    # Define parameters
    numberOfDimensions = 2
    numBins = len(data)

    r = np.linspace(0, maxRadius[i], numBins)

    r_mat = np.zeros((numBins, numBins))
    for j in range(numBins):
        for k in range(numBins):
            r_mat[j,k] = (j*j + k*k)**((numberOfDimensions-1)/2)

    data /= r_mat
    data /= np.sum(np.nan_to_num(data))

    plt.figure()
    plt.imshow(data, cmap=plt.cm.jet, extent=[0,maxRadius[i],maxRadius[i],0])
    plt.colorbar()
    plt.xlabel("Radial distance from particle 1")
    plt.ylabel("Radial distance from particle 2")
plt.show()
