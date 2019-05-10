import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# Load file
fileName = "../data/int1/twobody/RBMPJ/2D/12P/0.100000w/ADAM_MC262144.dat"
data = np.loadtxt(fileName)

# Define parameters
numberOfDimensions = 2
maxRadius = 25
numBins = len(data)

r = np.linspace(0, maxRadius, numBins)

r_mat = np.zeros((numBins, numBins))
for i in range(numBins):
    for j in range(numBins):
        r_mat[i,j] = (i*i + j*j)**((numberOfDimensions-1)/2)

data /= r_mat
data /= np.sum(np.nan_to_num(data))

plt.imshow(data, cmap=plt.cm.jet, extent=[0,maxRadius,maxRadius,0])
plt.xlabel("Radial distance from particle 1")
plt.ylabel("Radial distance from particle 2")
plt.show()
