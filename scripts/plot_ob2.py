import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm 
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from savitzky_golay import sgolay2d

plt.style.use("bmh")

class PlotOB():

    def __init__(self, fileName, radius):
        self.data = np.loadtxt(fileName)
        self.radius = radius
        
    def remove_cross(self):
        lengthHalf = int(len(self.data)/2)
        self.data[lengthHalf] = self.data[lengthHalf-1]
        self.data[lengthHalf+1] = self.data[lengthHalf+2]
        self.data[:,lengthHalf] = self.data[:,lengthHalf-1]
        self.data[:,lengthHalf+1] = self.data[:,lengthHalf+2]
    
    def norm(self, numberOfParticles):
        self.data = numberOfParticles * self.data/self.data.sum()
        self.data *= 1e4
        
    def cut(self, threshold=0.00001):
        self.data = np.where(self.data > threshold, 0, self.data)
        
    def crop(self, newRadius):
        length = len(self.data)
        newLength = length * (newRadius / self.radius)
        
        start = int(length/2-newLength/2)
        stop = int(length/2+newLength/2)
        
        self.data = self.data[start:stop, start:stop]
        self.radius = newRadius
        
    def rotate(self):
        self.data = np.rot90(self.data)
        
    def smooth(self, window, order):
        self.data = sgolay2d(self.data, window, order)
        
    def mask_cross_section(self):
        lengthHalf = int(len(self.data)/2)
        masked_data = self.data.copy()
        masked_data[:lengthHalf,lengthHalf:] = 0
        return masked_data
        
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
                    radiii[radii] += self.data[i, j]
                    
        plt.figure()
        x = np.linspace(0,self.radius,numberOfPoints/2)
        plt.plot(x, np.divide(radiii, binshit))
        plt.show()
        
    def fmt(self, x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
        
    def plot_3Dcontour(self, masked_data=None):
        if type(masked_data) is not np.ndarray:
            masked_data = self.data
        size = 24
        size_ticks = 24
        label_size = {"size":str(size)}
        plt.rcParams["font.family"] = "Serif"

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_facecolor('white')
        
        maximum = np.max(self.data)
        x = y = np.linspace(-self.radius, self.radius, len(self.data))
        xx, yy = np.meshgrid(x, y)

        ax.plot_surface(xx, yy, masked_data, rstride=8, cstride=8, alpha=1, antialiased=False)

        ax.contourf(xx, yy, self.data, zdir='z', offset=-maximum, cmap=cm.Blues)
        ax.contour(xx, yy, self.data, zdir='x', offset=-self.radius, levels=[0], cmap=cm.gist_gray)
        #ax.contourf(xx, yy, self.data, zdir='y', offset=self.radius)#, cmap=cm.gist_gray)
        

        ax.set_xlim(-self.radius, self.radius)
        ax.set_ylim(-self.radius, self.radius)
        ax.set_zlim(-maximum, maximum)

        ax.set_xlabel('$x$', labelpad=10, **label_size)
        ax.set_ylabel('$y$', labelpad=10, **label_size)
        ax.set_zlabel(r'$\rho(x,y)$', labelpad=20, **label_size)
        
        ticks = [-3, -2, -1, 0, 1, 2, 3]
        
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        #ax.zaxis.set_major_formatter(ticker.FuncFormatter(self.fmt))
        ax.zaxis.set_tick_params(pad=10)
        #ax.tick_params(labelsize=size_ticks)
        
        plt.tight_layout(rect=[-0.1, 0, 0.95, 1.06])
        plt.savefig("../plots/int0/onebody2/2D/2P/1.000000w/VMC_GD_MC2pow30_blue_small.png")
        plt.show()

if __name__ == "__main__":

    QD = PlotOB("../data/int0/quantumdot/onebody2/VMC/2D/2P/1.000000w/GD_MC65536.dat", 3)
    QD.remove_cross()
    QD.norm(6)
    #QD.cut(1)
    #QD.crop(50)
    #QD.smooth(29, 4)
    #masked_data = QD.mask_cross_section()
    #QD.plot_heatmap()
    #QD.rotate()
    #QD.plot_radial()
    QD.plot_3Dcontour()
