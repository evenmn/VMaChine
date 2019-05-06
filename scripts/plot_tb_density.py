import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# Load file
fileName = "../data/int1/twobody/VMC/2D/2P/1.000000w/SGD_MC262144_0.dat"
data = np.loadtxt(fileName)

# Define parameters
numberOfDimensions = 2
maxRadius = 5
numBins = len(data)

r = np.linspace(0, maxRadius, numBins)

r_mat = np.zeros((numBins, numBins))
for i in range(numBins):
    for j in range(numBins):
        r_mat[i,j] = (i*i + j*j)**((numberOfDimensions-1)/2)

data /= r_mat
data /= np.sum(np.nan_to_num(data))

plt.imshow(data, cmap=plt.cm.jet, extent=[0,5,5,0])
plt.show()
