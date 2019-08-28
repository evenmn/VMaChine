import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def exact(r1, r2, w):
    return np.exp(- w * (r1 * r1 + r2 * r2))
    
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
if __name__ == "__main__":
    N = 1000
    radius = 3
    r = np.linspace(-radius, radius, N)
    
    data = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            data[i, j] = exact(r[i], r[j], 1)
    
    data /= np.sum(data)
    
    size = 24
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams.update({'figure.autolayout': True})

    fig, ax = plt.subplots(figsize=(8,6))
    
    img = ax.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    cbar = fig.colorbar(img, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt))
    cbar.set_label(r'$\rho(r_i,r_j)$', rotation=270, labelpad=40, y=0.45, **label_size)
    
    plt.tight_layout()
    
    ax.set_xlabel("$r_j$", **label_size)
    ax.set_ylabel("$r_i$", **label_size)
    plt.show()
