import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../data/test.dat")

data /= data.sum()

radius = 5

plt.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
plt.colorbar()
plt.show()

'''
def doublewell(x, y, b):
    return x*x+y*y+b*b/4-b*abs(x)
    

x = np.arange(-10, 10, 0.1)
y = np.arange(-10, 10, 0.1)
xx, yy = np.meshgrid(x, y)
z = doublewell(xx,yy,4)
h = plt.contourf(x,y,z)
plt.show()
'''
