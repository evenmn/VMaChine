import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm 
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from savitzky_golay import savitzky_golay, sgolay2d

plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"
plt.rc('figure', max_open_warning = 0)
ax = plt.gca()
ax.set_facecolor('white')

class PlotED():
    def __init__(self, system = 'quantumdot', 
                       observable = 'onebody2',
                       method = 'VMC', 
                       dimension = 2, 
                       particles = 2, 
                       omega = '1.000000',
                       radius = 10,
                       interaction=1, 
                       optimizer='ADAM', 
                       cycles=1048576):
        '''
        Plotting engine for electron densities
        
        Arguments:
        ----------
        
        system      {string} :   System/potential
        observable  {string} :   Which electron density to plot
        methods     {string} :   Method
        dimensions  {int} :      Number of dimensions
        particles   {int} :      Number of particles
        omega       {string} :   Frequency
        radius      {float} :    Max radius of simulation
        newRadius   {float} :    New cropping radius
        interaction {bool} :     Interaction on (1)/ off (0)
        optimizer   {string} :   Optimizer used
        cycles      {int} :      Number of cycles, integer
        '''
        
        self.system = system
        self.method = method
        self.dimension = dimension
        self.particles = particles
        self.omega = omega
        self.radius = radius
        self.interaction = interaction
        self.optimizer = optimizer
        self.cycles = cycles
        
        fileName = self.generateFileName()
        self.load(fileName)
        
    def generateFileName(self):
        fileName  = "../data/int"
        fileName += str(self.interaction) + "/"
        fileName += self.system   + "/"
        fileName += "onebody2" + "/"
        fileName += self.method   + "/"
        fileName += str(self.dimension)      + "D/"
        fileName += str(self.particles) + "P/"
        fileName += str(self.omega)    + "w/"
        fileName += self.optimizer + "_MC"
        fileName += str(self.cycles) + ".dat"
        return fileName


    def saveFigure(self):
        fileName  = "../plots/int"
        fileName += str(self.interaction) + "/"
        fileName += "onebody2" + "/"
        fileName += str(self.dimension)      + "D/"
        fileName += str(self.particles) + "P/"
        fileName += str(self.omega)    + "w/"
        fileName += self.method + "_"  
        fileName += self.optimizer + "_MC"
        fileName += str(self.cycles) + ".png"
        print(fileName)
        plt.savefig(fileName)
        
    def load(self, fileName):
        self.data = np.loadtxt(fileName)
        
    def remove_cross(self):
        lengthHalf = int(len(self.data)/2)
        self.data[lengthHalf] = self.data[lengthHalf-1]
        self.data[lengthHalf+1] = self.data[lengthHalf+2]
        self.data[:,lengthHalf] = self.data[:,lengthHalf-1]
        self.data[:,lengthHalf+1] = self.data[:,lengthHalf+2]
        
    def radial(self):
        n = len(self.data)
        if self.dimension == 2:
            return 2*n+1
        elif self.dimension == 3:
            return 3*n*(n+1)+1
        else:
            raise ValueError("Number of dimensions needs to be 2 or 3!!!")
    
    def norm_radial(self, factor=1e4):
        v = self.radial()
        if len(self.data[0]) > 1:
            xx, yy = np.meshgrid(v, v, sparse=True)
            self.data /= np.multiply(xx,yy)
        else:
            self.data /= v
        self.data /= np.sum(np.nan_to_num(self.data))
        self.data *= factor * self.particles
    
    def norm(self, factor=1e4):
        norm_const = np.trapz(self.data)
        self.data /= norm_const
        self.data /= np.sum(np.nan_to_num(self.data))
        self.data *= factor * self.particles
        
        
    def cut(self, threshold=0.00001):
        self.data = np.where(self.data > threshold, 0, self.data)
        
    def crop_center(self, newRadius):
        length = len(self.data)
        newLength = length * (newRadius / self.radius)
        start = int(length/2-newLength/2)
        stop = int(length/2+newLength/2)
        self.data = self.data[start:stop, start:stop]
        self.radius = newRadius
        
    def crop_edges(self, newRadius):
        length = len(self.data)
        newLength = int(length * (newRadius / self.radius))
        self.data = self.data[0:newLength, 0:newLength]
        
    def window(self):
        self.data = self.data[2:,2:]

        data1 = np.rot90(self.data,0)
        data2 = np.rot90(self.data,-1)
        data3 = np.rot90(self.data,-2)
        data4 = np.rot90(self.data,-3)
        
        dataTop = np.c_[data2,data1]
        dataBottom = np.c_[data3,data4]
        self.data = np.concatenate((dataBottom,dataTop))
        
    def rotate(self):
        self.data = np.rot90(self.data)
        
    def smooth(self, window=29, order=4):
        if len(self.data[0]) > 1:
            self.data = sgolay2d(self.data, window, order)
        else:
            self.data = savitzky_golay(self.data, window, order)
        
    def mask_cross_section(self):
        lengthHalf = int(len(self.data)/2)
        masked_data = self.data.copy()
        masked_data[:lengthHalf,lengthHalf:] = 0
        return masked_data
        
    def ticks(self):
        if self.radius % 2 == 0 or self.radius < 2:
            return [-self.radius, -self.radius/2., 0, self.radius/2., self.radius]
        elif self.radius % 3 == 0:
            return [-self.radius, -2*(self.radius/3.), -self.radius/3., 0, self.radius/3., 2*(self.radius/3.), self.radius]
        elif self.radius == 5:
            return [-5, -3, -1, 1, 3, 5]
        elif self.radius == 25:
            return [-25, -15, -5, 5, 15, 25]
        else:
            print("Warning: Ticks out of bounds")
            return []
            
    def radial_profile(self):
        numberOfPoints = len(self.data)
        
        radiii = np.zeros(int(numberOfPoints/2))
        binshit = radiii.copy()
        for i in range(numberOfPoints):
            for j in range(numberOfPoints):
                radii = int(np.sqrt((i-numberOfPoints/2-1)**2 + (j-numberOfPoints/2-1)**2))
                if radii < numberOfPoints/2:
                    binshit[radii] += 1.
                    radiii[radii] += self.data[i, j]
                    
        x = np.linspace(0,self.radius,numberOfPoints/2)
        y = np.divide(radiii, binshit)
        return x, y
            
    def fmt(self, x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
        
    def plot_heatmap(self, ticks=None, 
                           size=24, 
                           size_ticks=14,
                           show=False,
                           save=False):
        '''
        Plot radial two-body density heatmap profile.
        
        Arguments:
        ---------
        
        ticks       {list} :        Axis ticks
        size        {int} :         Label size
        size_ticks  {int} :         Ticks size
        show        {bool} :        Display figure
        save        {bool} :        Save figure
        '''
        if type(ticks) is not list:
            ticks = self.ticks()
    
        label_size = {"size":str(size)}
        plt.rcParams["font.family"] = "Serif"
        plt.rcParams.update({'figure.autolayout': True})

        fig, ax = plt.subplots(figsize=(8,6))
                     
        img = ax.imshow(self.data, cmap=plt.cm.jet, extent=[-self.radius,self.radius,-self.radius,self.radius])
        cbar = fig.colorbar(img, fraction=0.046, pad=0.04)
        cbar.set_label(r'$\rho(r_i,r_j)$', rotation=270, labelpad=40, y=0.45, **label_size)
        cbar.ax.tick_params(labelsize=size_ticks)
        
        plt.tight_layout()
        
        ax.set_xlabel("$r_j$", **label_size)
        ax.set_ylabel("$r_i$", **label_size)
        ax.tick_params(labelsize=size_ticks)
            
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        plt.grid()
        
        if save:
            self.saveFigure()
        if show:
            plt.show()  
        
    def plot_radial(self, label,
                          show=False,
                          save=False):
        '''
        Plot radial one-body density profile.
        
        Arguments:
        ---------
        
        show        {bool} :        Display figure
        save        {bool} :        Save figure
        '''
        
        x, y = self.radial_profile()
        plt.plot(x, y, label=label)
        if save:
            self.saveFigure()
        if show:
            plt.show()  
        
    def plot_3D(self, ticks=None, 
                      masked_data=None, 
                      size=24, 
                      size_ticks=14,
                      show=False,
                      save=False):
        '''
        Plot 3D surface.
        
        Arguments:
        ---------
        
        ticks       {list} :        Axis ticks
        masked_data {np.ndarray} :  Masked data to be used in surface plot
        size        {int} :         Label size
        size_ticks  {int} :         Ticks size
        show        {bool} :        Display figure
        save        {bool} :        Save figure
        '''
                      
        if type(masked_data) is not np.ndarray:
            masked_data = self.data
        if type(ticks) is not list:
            ticks = self.ticks()
    
        label_size = {"size":str(size)}

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_facecolor('white')
        
        maximum = np.max(self.data)
        v = np.linspace(-self.radius, self.radius, len(self.data))
        xx, yy = np.meshgrid(v, v)

        surf = ax.plot_surface(xx, yy, masked_data, rstride=8, cstride=8, alpha=1, antialiased=False)

        ax.set_xlim(-self.radius, self.radius)
        ax.set_ylim(-self.radius, self.radius)
        ax.set_zlim(0, maximum)
        
        ax.set_xlabel('$x$', labelpad=10, **label_size)
        ax.set_ylabel('$y$', labelpad=10, **label_size)
        ax.set_zlabel(r'$\rho(x,y)$', labelpad=20, **label_size)
        
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.zaxis.set_tick_params(pad=10)
        ax.tick_params(labelsize=size_ticks)
        
        plt.tight_layout(rect=[-0.1, 0, 0.95, 1.06])
        if save:
            self.saveFigure()
        if show:
            plt.show()  
        
    def plot_3Dcontour(self, ticks=None, 
                             masked_data=None, 
                             size=24, 
                             size_ticks=12,
                             show=False,
                             save=False):
        '''
        Plot 3D surface with contour plot on the xy-plane
        and cross section line on the yz-plane.
        
        Arguments:
        ---------
        
        ticks       {list} :        Axis ticks
        masked_data {np.ndarray} :  Masked data to be used in surface plot
        size        {int} :         Label size
        size_ticks  {int} :         Ticks size
        show        {bool} :        Display figure
        save        {bool} :        Save figure
        '''
                             
        if type(masked_data) is not np.ndarray:
            masked_data = self.data
        if type(ticks) is not list:
            ticks = self.ticks()
            
        label_size = {"size":str(size)}

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_facecolor('white')
        
        maximum = np.max(self.data)
        v = np.linspace(-self.radius, self.radius, len(self.data))
        xx, yy = np.meshgrid(v, v)

        ax.plot_surface(xx, yy, masked_data, rstride=8, cstride=8, alpha=1, antialiased=False)
        ax.contour(xx, yy, self.data, zdir='x', offset=-self.radius, levels=[0], cmap=cm.gist_gray)
        ax.contourf(xx, yy, self.data, zdir='z', offset=-maximum, cmap=cm.Blues)
        #ax.contourf(xx, yy, self.data, zdir='y', offset=self.radius)#, cmap=cm.gist_gray)
        

        ax.set_xlim(-self.radius, self.radius)
        ax.set_ylim(-self.radius, self.radius)
        ax.set_zlim(-maximum, maximum)

        ax.set_xlabel('$x$', labelpad=10, **label_size)
        ax.set_ylabel('$y$', labelpad=10, **label_size)
        ax.set_zlabel(r'$\rho(x,y)$', labelpad=20, **label_size)
        
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.zaxis.set_tick_params(pad=10)
        ax.tick_params(labelsize=size_ticks)
        
        plt.tight_layout(rect=[-0.1, 0, 0.95, 1.06])
        if save:
            self.saveFigure()
        if show:
            plt.show()     
        
def saveFigure(dim, par, omega, optimizer='ADAM', cycles=1048576):
        fileName  = "../plots/int1/"
        fileName += "onebody" + "/"
        fileName += str(dim)      + "D/"
        fileName += str(par) + "P/"
        fileName += str(omega)    + "w/"
        fileName += optimizer + "_MC"
        fileName += str(cycles) + ".png"
        print(fileName)
        plt.savefig(fileName)

if __name__ == "__main__":

    system   = 'quantumdot'
    observable = 'onebody2'
    
    methods   = [#'VMC',
                 'RBM',
                 #'RBMSJ',
                 #'RBMPJ'
                ]
                
    dims      = [2]
    
    particles = [2, 
                 #6, 
                 #12, 
                 #20
                 #30,
                 #42,
                 #56,
                 #2,
                 #8,
                 #20
                 ]
                 
    omegas    = ['1.000000',
                 #'0.500000',
                 #'0.280000',
                 #'0.100000',
                 #'0.010000'
                 ] 

    radius = [10,
              #15, 
              #25, 
              #30,
              #35,
              #40,
              #45
              #50
              ]
                 
    newRadius = [#3, 4, 6, 10,
                 #4, 6, 8, 15,
                 #5, 8, 10, 20,
                 #6, 10, 12, 25,
                 #6, 12, 14, 30,
                 #7, 14, 16, 35
                 3
                 ]
                 
    for d in range(len(dims)):
        for p in range(len(particles)):
            for o in range(len(omegas)):
                for m in range(len(methods)):
                    QD = PlotED(system, observable, methods[m], dims[d], particles[p], omegas[o], radius[p], interaction=0, cycles=268435456)
                    QD.crop_center(newRadius[p])
                    QD.norm_radial()
                    QD.remove_cross()
                    #QD.cut(0.2)
                    QD.smooth()
                    QD.plot_3Dcontour(save=True, show=True)
    
