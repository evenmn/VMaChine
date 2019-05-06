import numpy as np
import matplotlib.pyplot as plt

infile = "../data/int1/onebody/VMC/2D/2P/1.000000w/SGD_MC262144_0.dat"
data = np.loadtxt(infile)

print(np.sum(data))
