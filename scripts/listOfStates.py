import numpy as np

orbitals = 5

for i in range(orbitals):
    for j in range(i+1):
        print(i-j,j)
