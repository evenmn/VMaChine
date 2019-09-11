import numpy as np
import matplotlib.pyplot as plt
from savitzky_golay import savitzky_golay

plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"
plt.rc('figure', max_open_warning = 0)

class PlotOB():
    def __init__(self, system = 'quantumdot', 
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
        methods     {string} :   Method
        dimensions  {int} :      Number of dimensions
        particles   {int} :      Number of particles
        omega       {float} :    Frequency
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
        self.r = np.linspace(0, self.radius, len(self.data))
        
    def generateFileName(self):
        fileName  = "../data/int"
        fileName += str(self.interaction) + "/"
        fileName += self.system   + "/"
        fileName += "onebody" + "/"
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
        fileName += "onebody" + "/"
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
        
    def radial(self):
        n = np.arange(len(self.data))
        if self.dimension == 2:
            return 2*n+1
        elif self.dimension == 3:
            return 3*n*(n+1)+1
        else:
            return np.nan
        
    def norm_bins(self):
        self.data /= self.radial()
        
    def norm(self):
        #data /= np.sum(np.nan_to_num(data))
        #data /= maxRadius[i]
        #data /= np.sum(data)
        #data *= int(len(data)/1000)
        self.data /= self.cycles
    
        #data = np.multiply(data, r**2)
        #r = r[np.argmax(data):]
        #data = data[np.argmax(data):]
        norm_const = np.trapz(self.data)
        self.data /= norm_const
        self.data *= self.particles
        self.data *= 100
        if self.method == 'VMC':
            self.data/=3
        
    def crop_edges(self, newRadius):
        length = len(self.data)
        newLength = int(length * (newRadius / self.radius))
        self.data = self.data[0:newLength]
        self.radius = newRadius
        self.r = np.linspace(0, self.radius, len(self.data))
        
    def cut(self, threshold=0.00001):
        self.data = np.where(self.data > threshold, 0, self.data)
        
    def cut_noise(self, limit):
        self.data[:np.argmax(self.data)] = np.where(self.data[:np.argmax(self.data)] < limit, 0, self.data[:np.argmax(self.data)])
        indices = np.where(self.data == 0)[0]
        self.data = np.delete(self.data, indices)
        self.r = np.delete(self.r, indices)
        
    def smooth(self, window=29, order=4):
        if len(self.data[0]) > 1:
            self.data = sgolay2d(self.data, window, order)
        else:
            self.data = savitzky_golay(self.data, window, order)

    def exact(self):
        '''Exact solution without interaction for given w'''
        return np.sqrt(self.omega/np.pi) * np.exp(- self.omega * self.r**2)
            
    def plot_radial(self, line_style, label):
        plt.plot(self.r, self.data, line_style, markersize=1, label=label)
        
        
def saveFigure(dim, par, omega, method, optimizer='ADAM', cycles=1048576):
        fileName  = "../plots/int1/"
        fileName += "onebody" + "/"
        fileName += str(dim)      + "D/"
        fileName += str(par) + "P/"
        fileName += str(omega)    + "w/"
        fileName += optimizer + "_MC"
        fileName += str(cycles) + ".png"
        print(fileName)
        plt.savefig(fileName)
        
if __name__ == '__main__':

    system   = 'quantumdot'
    
    methods   = ['VMC',
                 'RBM',
                 'RBMSJ',
                 #'RBMPJ'
                ]
                
    dims      = [2]
    
    particles = [42, 
                 #6, 
                 #12, 
                 #20
                 #30,
                 #56,
                 #2,
                 #8,
                 #20
                 ]
                 
    omegas    = [#'1.000000',
                 #'0.500000',
                 #'0.280000',
                 '0.100000',
                 #'0.010000'
                 ] 

    radius = [#10, 
              #15, 
              #25,
              #30,
              #35,
              40,
              #55
              ]
                 
    newRadius = [#3, 4, 4, 10,
                 #4, 6, 8, 15,
                 #5, 8, 10, 20,
                 #6, 8, 12, 25,
                 #7, 8, 14, 25,
                 #7, 10, 16, 25,
                 25
                 ]
                 
    line_style = ["-", 
                  "--", 
                  "-.", 
                  ":"
                  ]
                  
    limits = [0.0, 0.0, 0.0, 0.0]
           
    i = 0
    for d in range(len(dims)):
        for p in range(len(particles)):
            for o in range(len(omegas)):
                plt.figure()
                ax = plt.gca()
                ax.set_facecolor('white')
                ax.tick_params(labelsize=16)
                for m in range(len(methods)):
                    QD = PlotOB(system, methods[m], dims[d], particles[p], omegas[o], radius[p])
                    QD.crop_edges(newRadius[i])
                    QD.norm_bins()
                    QD.norm()
                    QD.cut_noise(limits[m])
                    #QD.smooth()
                    QD.plot_radial(line_style=line_style[m], label=methods[m])
                  
                size = 24
                size_legend=16
                label_size = {"size":str(size)}

                plt.gcf().subplots_adjust(bottom=0.15)
                plt.gcf().subplots_adjust(left=0.18)
                #plt.tight_layout()

                plt.xlabel("$r$", **label_size)
                plt.ylabel(r"$\rho(r)$", **label_size)
                plt.legend(loc="best", fontsize=size_legend, facecolor='white', framealpha=1)
                saveFigure(dims[d], particles[p], omegas[o], methods[m])
                plt.show()
                i += 1
