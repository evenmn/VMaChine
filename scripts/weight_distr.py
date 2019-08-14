import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../data/int1/weights/RBMPJ/2D/6P/1.000000w/ADAM_MC1048576.dat')

data_prod = data[2,:]
data_gauss = data[1,:]
data_a = data_gauss[:6]
data_b = data_prod[:6]

data_a = data_a.reshape((1,-1))
data_b = data_b.reshape((1,-1))

data_w = data_prod[-72:]
data_w_reshaped = data_w.reshape((12,6))

print(data_w)
print(data_w_reshaped)

fig, ax = plt.subplots(ncols=3)
im1 = ax[0].imshow(data_a.T)
im2 = ax[1].imshow(data_b.T)
im3 = ax[2].imshow(data_w_reshaped)
#plt.colorbar()
plt.show()
