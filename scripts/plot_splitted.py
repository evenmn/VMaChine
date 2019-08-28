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
    size = 16
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.18)

    plt.plot(w, A, label="VMC")
    plt.plot(w, B, label="RBM")
    plt.plot(w, C, label="RBM+SJ")
    plt.plot(w, D, label="RBM+PJ")
    plt.xlabel("$\omega$", **label_size)
    plt.ylabel(r"$\frac{\langle\hat{T}\rangle}{\langle\hat{H}\rangle}$", rotation=0, labelpad=10, **label_size)
    plt.legend(loc="best", fontsize=size)
    plt.grid()
    plt.show()
    
if __name__ == "__main__":
    w = [0.01,0.1,0.28,0.5,1.0,2.0,3.0,5.0,10.0]
    RBMPJ_2D_2P = [[0.074107, 0.440975, 1.021668, 1.659637, 2.999587, 5.49475, 7.87961, 12.49832, 23.65070],
                   [0.01031, 0.09223, 0.2468, 0.4305, 0.8440, 1.7234, 2.3144, 3.9569, 5.384],
                   [0.02703, 0.1757, 0.4258, 0.7112, 1.3418, 2.4657, 3.9349, 6.3068, 15.265],
                   [0.03677, 0.17304, 0.3490, 0.5179, 0.8238, 1.3057, 1.6413, 2.2347, 3.0010]]
                   
    VMC_2D_2P = [[0.074070, 0.44129, 1.02192, 1.65974, 2.99936, 5.497567, 7.881906, 12.504211, 23.65035],
                 [0.00947, 0.09117, 0.2477, 0.4346, 0.8523, 2.0574, 2.5028, 5.1819, 9.529],
                 [0.02732, 0.1789, 0.4256, 0.7057, 1.3149, 1.9931, 3.6911, 4.7651, 10.538],
                 [0.03728, 0.17119, 0.3487, 0.5195, 0.8321, 1.4471, 1.6880, 2.5572, 3.583]]
                 
    RBMSJ_2D_2P = [[0.075267, 0.44858, 1.03470, 1.67636, 3.0213, 5.5254, 7.9103, 12.5322, 23.7539],
                   [0.00738, 0.07539, 0.2163, 0.3967, 0.8120, 1.6572, 2.5627, 4.369, 5.455],
                   [0.03071, 0.1990, 0.4547, 0.7366, 1.3623, 2.5447, 3.664, 5.890, 15.279],
                   [0.03718, 0.1742, 0.3637, 0.5430, 0.8469, 1.3234, 1.6840, 2.2838, 3.0193]]
                   
    RBM_2D_2P = [[0.07954, 0.4743, 1.0707, 1.7234, 3.0829, 5.5936, 7.9968, 12.6070, 23.6589],
                 [0.00872, 0.08102, 0.2047, 0.3739, 0.7691, 1.6377, 2.2346, 3.8768, 5.054],
                 [0.03402, 0.2082, 0.4678, 0.7611, 1.3926, 2.5507, 4.0201, 6.4292, 15.592],
                 [0.0368, 0.1851, 0.3983, 0.5884, 0.9212, 1.4051, 1.7422, 2.3010, 3.0129]]
            
    #plot(w, RBMPJ_2D_2P)       
    #plot(w, VMC_2D_2P)
    
    A = np.divide(VMC_2D_2P[3],VMC_2D_2P[0]) 
    B = np.divide(RBM_2D_2P[3],RBMPJ_2D_2P[0]) 
    C = np.divide(RBMSJ_2D_2P[3],RBMPJ_2D_2P[0]) 
    D = np.divide(RBMPJ_2D_2P[3],RBMPJ_2D_2P[0]) 
    
    plot_split(A, B, C, D, w)
    
    #tikzplotlib.save("/home/evenmn/Master-thesis/doc/text/pgf/kinetic.tex")
