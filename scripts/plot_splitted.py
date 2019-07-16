import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
            
def plot(w, data):  
    kin = np.divide(data[1],data[0]) 
    ext = np.divide(data[2],data[0]) 
    inr = np.divide(data[3],data[0]) 
    plt.plot(w, kin, label=r"$\langle\mathcal{T}\rangle/\langle\mathcal{H}\rangle$")
    plt.plot(w, ext, label=r"$\langle\mathcal{V}_{ext}\rangle/\langle\mathcal{H}\rangle$")
    plt.plot(w, inr, label=r"$\langle\mathcal{V}_{int}\rangle/\langle\mathcal{H}\rangle$")
    plt.xlabel("$\omega$")
    plt.ylabel("Energy")
    plt.legend(loc="best")
    plt.grid()
    #plt.show()
    
def plot_split(A, B, C, D, w):
    plt.plot(w, A, label="VMC")
    plt.plot(w, B, label="RBM")
    plt.plot(w, C, label="RBM+SJ")
    plt.plot(w, D, label="RBM+PJ")
    plt.xlabel("$\omega$")
    plt.ylabel(r"$\langle\mathcal{T}\rangle/\langle\mathcal{H}\rangle$")
    plt.legend(loc="best")
    plt.grid()
    #plt.show()
    
if __name__ == "__main__":
    w = [0.01,0.1,0.28,0.5,1.0,3.0]
    RBMPJ_2D_2P = [[0.074107, 0.440975, 1.021668, 1.659637, 2.999587, 7.882355],
                   [0.01031, 0.09223, 0.2468, 0.4305, 0.8440, 2.2466],
                   [0.02703, 0.1757, 0.4258, 0.7112, 1.3418, 4.0117],
                   [0.03677, 0.17304, 0.3490, 0.5179, 0.8238, 1.6240]]
                   
    VMC_2D_2P = [[0.074070, 0.44129, 1.02192, 1.65974, 2.99936, 7.881906],
                 [0.00947, 0.09117, 0.2477, 0.4346, 0.8523, 2.5028],
                 [0.02732, 0.1789, 0.4256, 0.7057, 1.3149, 3.6911],
                 [0.03728, 0.17119, 0.3487, 0.5195, 0.8321, 1.6880]]
                 
    RBMSJ_2D_2P = [[0.075267, 0.44858, 1.03470, 1.167739, 3.0259, 7.9409],
                   [0.00738, 0.07539, 0.2163, 0.3913, 0.7857, 2.2162],
                   [0.03071, 0.1990, 0.4547, 0.7450, 1.3958, 4.0834],
                   [0.03718, 0.1742, 0.3637, 0.5411, 0.8444, 1.6413]]
                   
    RBM_2D_2P = [[0.07954, 0.4743, 1.0707, 1.7234, 3.0829, 7.9968],
                 [0.00872, 0.08102, 0.2047, 0.3739, 0.7691, 2.2346],
                 [0.03402, 0.2082, 0.4678, 0.7611, 1.3926, 4.0201],
                 [0.0368, 0.1851, 0.3983, 0.5884, 0.9212, 1.7422]]
            
    #plot(w, RBMPJ_2D_2P)       
    #plot(w, VMC_2D_2P)
    
    A = np.divide(VMC_2D_2P[1],VMC_2D_2P[0]) 
    B = np.divide(RBM_2D_2P[1],RBMPJ_2D_2P[0]) 
    C = np.divide(RBMSJ_2D_2P[1],RBMPJ_2D_2P[0]) 
    D = np.divide(RBMPJ_2D_2P[1],RBMPJ_2D_2P[0]) 
    
    plot_split(A, B, C, D, w)
    
    tikzplotlib.save("/home/evenmn/Master-thesis/doc/text/pgf/kinetic.tex")
