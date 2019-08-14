import numpy as np

N = 3
D = 3

G = np.zeros((N*D,N*D))

for n in range(1, N):
    for j in range(D*(N-n)):
        G[j, j+n*D] = j + 10 * (j+n*D)
        
print(G)
