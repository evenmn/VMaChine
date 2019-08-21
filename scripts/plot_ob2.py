import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm 
from matplotlib.ticker import LinearLocator, FormatStrFormatter

class PlotOB():

    def __init__(self, fileName, radius):
        self.data = np.loadtxt(fileName)
        self.radius = radius
    
    def norm(self, numberOfParticles):
        self.data = numberOfParticles * self.data/self.data.sum()
        
    def cut(self, threshold=0.00001):
        self.data = np.where(self.data > threshold, 0, self.data)
        
    def crop(self, newRadius):
        length = len(self.data)
        newLength = length * (newRadius / self.radius)
        
        start = int(length/2-newLength/2)
        stop = int(length/2+newLength/2)
        
        self.data = self.data[start:stop, start:stop]
        self.radius = newRadius
        
    def plot_heatmap(self):
        plt.figure()
        plt.imshow(self.data, cmap=plt.cm.jet, extent=[-self.radius,self.radius,
                                                       -self.radius,self.radius])
        plt.colorbar()
        plt.show()
        
    def plot_3D(self):
        x = y = np.linspace(-self.radius, self.radius, len(self.data))
        xx, yy = np.meshgrid(x, y, sparse=True)

        fig = plt.figure() 
        ax = fig.gca(projection='3d')

        #Plot the surface. 
        surf = ax.plot_surface(xx, yy, self.data, cmap=cm.coolwarm, linewidth=0, antialiased=False)

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
        
    def plot_radial(self):
        numberOfPoints = len(self.data)
        
        radiii = np.zeros(int(numberOfPoints/2))
        binshit = radiii.copy()
        for i in range(numberOfPoints):
            for j in range(numberOfPoints):
                radii = int(np.sqrt((i-numberOfPoints/2-1)**2 + (j-numberOfPoints/2-1)**2))
                if radii < numberOfPoints/2:
                    binshit[radii] += 1.
                    radiii[radii] += data[i, j]
                    
        plt.figure()
        x = np.linspace(0,self.radius,numberOfPoints/2)
        plt.plot(x, np.divide(radiii, binshit))
        plt.show()
        
    def fmt(self, x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
        
    def plot_3Dcontour(self):
        size = 16
        label_size = {"size":str(size)}
        plt.rcParams["font.family"] = "Serif"
        plt.rcParams['mathtext.default'] = 'regular'
        #plt.rcParams.update({'figure.autolayout': True})

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        maximum = np.max(self.data)
        x = y = np.linspace(-self.radius, self.radius, len(self.data))
        xx, yy = np.meshgrid(x, y)

        ax.plot_surface(xx, yy, self.data, rstride=8, cstride=8, alpha=0.5)

        cset = ax.contourf(xx, yy, self.data, zdir='z', offset=-maximum, cmap=cm.jet)
        cset = ax.contourf(xx, yy, self.data, zdir='x', offset=-self.radius, cmap=cm.coolwarm)
        cset = ax.contourf(xx, yy, self.data, zdir='y', offset=self.radius, cmap=cm.coolwarm)

        ax.set_xlim(-self.radius, self.radius)
        ax.set_ylim(-self.radius, self.radius)
        ax.set_zlim(-maximum, maximum)

        ax.set_xlabel('$x$', labelpad=10, **label_size)
        ax.set_ylabel('$y$', labelpad=10, **label_size)
        ax.set_zlabel(r'$\rho$(x,y)', labelpad=30, **label_size)
        
        ax.set_xticks([-3,-2,-1, 0, 1, 2, 3])
        ax.set_yticks([-3,-2,-1, 0, 1, 2, 3])
        ax.zaxis.set_major_formatter(ticker.FuncFormatter(self.fmt))
        ax.zaxis.set_tick_params(pad=15)
        
        plt.tight_layout(rect=[-0.1, 0, 0.91, 1.06])
        plt.show()

if __name__ == "__main__":

    QD = PlotOB("../data/int1/quantumdot/onebody2/VMC/2D/2P/1.000000w/ADAM_MC262144.dat", 5)
    QD.norm(2)
    QD.cut(0.000005)
    QD.crop(3)
    QD.plot_3Dcontour()
