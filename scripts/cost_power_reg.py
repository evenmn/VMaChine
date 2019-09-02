import numpy as np
import matplotlib.pyplot as plt

def reg_power(x, y):
    ''' Returns the coefficients a, b in f(x)=ax^b
    
    Based on regression formulas from
    http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
    '''
    
    n = len(x)
    lnx = 0
    lny = 0
    lnxlnx = 0
    lnxlny = 0
    for i in range(n):
        lnx += np.log(x[i])
        lny += np.log(y[i])
        lnxlnx += np.log(x[i])*np.log(x[i])
        lnxlny += np.log(x[i])*np.log(y[i])
    
    b = (n * lnxlny - lnx*lny) / (n * lnxlnx - lnx*lnx)
    a = (lny - b * lnx) / n
    return np.exp(a), b
    
def f(x, a, b):
    return a*x**b

if __name__ == '__main__':
    P_2D = [2, 6, 12, 20, 30, 42, 56, 72, 90]
    RBM_2D = [6.05, 11.25, 20.53, 38.99, 73.72, 130.49, 213.47, 360.22,856.84]
    RBMSJ_2D = [7.12, 14.07, 28.42, 63.27, 122.93, 199.60, 349.22]
    RBMPJ_2D = [7.26, 13.50, 27.68, 57.09, 119.17, 212.53, 382.13]
    VMC_2D = [5.11, 10.51, 20.85, 41.20, 76.26, 137.39, 230.63, 355.81, 544.03]

    P_3D = [2, 8, 20, 40, 70]
    RBM_3D = [7.69, 20.92, 59.67, 171.84, 586.39]
    RBMSJ_3D = [8.95, 26.86, 94.64, 270.92]
    RBMPJ_3D = [8.87, 26.36, 91.40, 293.25]
    VMC_3D = [6.70, 20.99, 62.54, 185.65, 486.02]
    
    a, b = reg_power(P_2D, VMC_2D)
    

    plt.plot(P_2D, VMC_2D)
    plt.plot(P_2D, f(P_2D, a, b))
    #plt.plot(P_2D[:-2], RBMSJ_2D)
    #plt.plot(P_2D[:-2], RBMPJ_2D)
    #plt.plot(P_2D, VMC_2D)
    plt.show()
    
    
