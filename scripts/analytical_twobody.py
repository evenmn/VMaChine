import numpy as np
import matplotlib.pyplot as plt

w = 1

def phi(r1, r2):
    return np.exp(-w*(r1*r1+r2*r2))

N = 1000

x = np.linspace(-2,2,N)

data = np.zeros((N, N))

i = 0
for r1 in x:
    j = 0
    for r2 in x:
        data[i,j] = phi(r1, r2)
        j += 1
    i += 1
    
plt.imshow(data)
plt.colorbar()
plt.show()
