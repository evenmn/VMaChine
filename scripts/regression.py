import numpy as np
import matplotlib.pyplot as plt

x_max = 5
N = 100

x = np.linspace(-x_max, x_max, N)

XXXX, YYYY, ZZZZ, VVVV = np.meshgrid(x, x, x, x, indexing='ij')

X = np.stack((XXXX, YYYY, ZZZZ, VVVV), axis=4)
print(X.shape)

def f(x):
    return np.exp(-0.5 * np.square(x))
