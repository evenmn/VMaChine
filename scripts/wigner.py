import numpy as np
import matplotlib.pyplot as plt

class Crystalization:
    def __init__(self, N, D, cycles, eta):
        self.N = N
        self.D = D
        self.cycles = cycles
        self.eta = eta
        self.initialize()
        
    def initialize(self):
        self.positions = np.random.rand(self.N, self.D)
        self.matrix = np.zeros((self.N, self.N))
        self.vector = np.zeros(self.N)
        self.radialVector()
        self.distanceMatrix()
        self.ext = self.vector.sum()/2
        self.int = self.matrix.sum()/2
        self.pot = self.ext + self.int

    def distanceMatrix(self):
        for i in range(self.N):
            for j in range(self.N):
                counter = 0
                for d in range(self.D):
                    counter += (self.positions[i, d] - self.positions[j, d])**2
                self.matrix[i, j] = counter
        
    def radialVector(self):
        for i in range(self.N):
            counter = 0
            for d in range(self.D):
                counter += self.positions[i, d]**2
            self.vector[i] = counter
            
    def set(self):
        self.oldPositions = self.positions
        self.oldVector = self.vector
        self.oldMatrix = self.matrix
        
    def reset(self):
        self.positions = self.oldPositions
        self.vector = self.oldVector
        self.matrix = self.oldMatrix

    def __iter__(self):
        for iter in range(self.cycles):
            self.set()
            self.positions += self.eta * np.random.rand(self.N, self.D)
            self.radialVector()
            self.distanceMatrix()
            self.ext = self.vector.sum()/2
            self.int = self.matrix.sum()/2
            #self.pot = self.ext + self.int
            
            if self.ext + self.int > self.pot:
                self.reset()
            else:
                self.pot = self.ext + self.int
                
    def plot(self):
        for i in range(self.N):
            plt.plot(self.positions[i,0], self.positions[i,1], 'or')
        plt.axis("equal")
        plt.show()
                
if __name__ == "__main__":
    N = 12
    D = 2
    cycles = 100000
    eta = 0.001
    
    crystal = Crystalization(N, D, cycles, eta)
    crystal.__iter__()
    crystal.plot()
                
