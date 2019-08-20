import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm 
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def prepare_data(fileName="../data/test.dat"):
    data = np.loadtxt(fileName)
    data /= data.sum()
    return np.where(data > 0.00002, 0, data)
    
def plot_heatmap(data, radius=5):
    plt.figure()
    plt.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    plt.colorbar()
    plt.show()
    
def plot_3D(data, radius=5):
    x = y = np.linspace(-radius, radius, 3000)
    xx, yy = np.meshgrid(x, y, sparse=True)

    fig = plt.figure() 
    ax = fig.gca(projection='3d')

    #Plot the surface. 
    surf = ax.plot_surface(xx,yy,data,cmap=cm.coolwarm,linewidth=0,antialiased=False)

    #Customize the z axis. 
    #ax.set_zlim(-0.10,1.40)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    #Add a color bar which maps values to colors. 
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    
def plot_radial(data, radius=5):
    numberOfPoints = len(data)
    radiii = np.zeros(int(numberOfPoints/2))
    binshit = radiii.copy()
    for i in range(numberOfPoints):
        for j in range(numberOfPoints):
            radii = int(np.sqrt((i-numberOfPoints/2-1)**2 + (j-numberOfPoints/2-1)**2))
            if radii < numberOfPoints/2:
                binshit[radii] += 1.
                radiii[radii] += data[i, j]
                
    plt.figure()
    x = np.linspace(0,radius,numberOfPoints/2)
    plt.plot(x, np.divide(radiii, binshit))
    plt.show()

if __name__ == "__main__":
    radius = 10

    data = prepare_data("../data/int1/quantumdot/onebody2/VMC/2D/2P/1.000000w/ADAM_MC262144.dat")
    plot_heatmap(data, radius)
    plot_radial(data, radius)
