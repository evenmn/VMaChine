import sympy as sp
import numpy as np

def factorial(n):
    if n <= 1:
        return 1
    else:
        return factorial(n - 1) * n

def symbolic_integration(n):
    x = sp.Symbol('x')
    ubestemt_integral = sp.integrate(x**n * sp.exp(-x*x))
    øvre  = ubestemt_integral.subs(x,100000)
    nedre = ubestemt_integral.subs(x,0)
    bestemt_integral = øvre - nedre
    return bestemt_integral
    
def iterative_integration(n):
    if(n%2 == 0):
        double_fact  = factorial(factorial(n-1))
        nominator = 2**(n/2 + 1)
        return (double_fact/nominator) * np.sqrt(np.pi)
    else:
        return factorial((n-1)/2)/2

def hermite_coefficients(n):
    if n == 0:
        return [1]
    elif n == 1:
        return [0,2]
    elif n == 2:
        return [-2,0,4]
    elif n == 3:
        return [0,-12,0,8]
    elif n == 4:
        return [12,0,-48,0,16]
    elif n == 5:
        return [0,120,0,-160,0,32]
    elif n == 6:
        return [-120,0,720,0,-480,0,64]
    elif n == 7:
        return [0,-1680,0,3360,0,-1344,0,128]
    elif n == 8:
        return [1680,0,-13440,0,13440,0,-3584,0,256]
    elif n == 9:
        return [0,30240,0,-80640,0,48384,0,-9326,0,512]
    elif n == 10:
        return [-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024]
    elif n == 11:
        return [0,-665280,0,2217600,0,-1774080,0,506880,0,-56320,0,2048]
    elif n == 12:
        return [665280,0,-7983360,0,13305600,0,-7096320,0,1520640,0,-135168,0,4096]
    else:
        return 0
        
'''
double H13(const double x) {return 8192*pow(x, 13) - 319488*pow(x, 11) + 4392960*pow(x, 9) - 26357760*pow(x, 7) + 69189120*pow(x, 5) - 69189120*pow(x, 3) + 17297280*x;}
double H14(const double x) {return 16384*pow(x, 14) - 745472*pow(x, 12) + 12300288*pow(x, 10) - 92252160*pow(x, 8) + 322882560*pow(x, 6) - 484323840*pow(x, 4) + 242161920*pow(x, 2) - 17297280;}
double H15(const double x) {return 32768*pow(x, 15) - 1720320*pow(x, 13) + 33546240*pow(x, 11) - 307507200*pow(x, 9) + 1383782400*pow(x, 7) - 2905943040*pow(x, 5) + 2421619200*pow(x, 3) - 518918400*x;}
double H16(const double x) {return 65536*pow(x, 16) - 3932160*pow(x, 14) + 89456640*pow(x, 12) - 984023040*pow(x, 10) + 5535129600*pow(x, 8) - 15498362880*pow(x, 6) + 19372953600*pow(x, 4) - 8302694400*pow(x, 2) + 518918400;}
double H17(const double x) {return 131072*pow(x, 17) - 8912896*pow(x, 15) + 233963520*pow(x, 13) - 3041525760*pow(x, 11) + 20910489600*pow(x, 9) - 75277762560*pow(x, 7) + 131736084480*pow(x, 5) - 94097203200*pow(x, 3) + 17643225600*x;}
'''

def merge_coefficient_lists(m,n):
    hermite_m = hermite_coefficients(m)
    hermite_n = hermite_coefficients(n)
    new_list = np.zeros(len(hermite_m) + len(hermite_n) - 1)
    for i in range(len(hermite_m)):
        for j in range(len(hermite_n)):
            new_list[i+j] += hermite_m[i] * hermite_n[j]
    return new_list

def solve_hermite_integrals(m,n):
    coeffs = merge_coefficient_lists(m,n)
    coeffs = np.insert(coeffs,0,0)
    summer = 0
    for i in range(len(coeffs)):
        summer += coeffs[i]*iterative_integration(i)
    return summer

def calculate_h_plus(basisSize, omega):
    h_plus = np.zeros((basisSize,basisSize))
    for m in range(basisSize):
        for n in range(basisSize):
            h_plus[m,n] = solve_hermite_integrals(m,n)
    return h_plus/omega
    
def calculate_h_HO(basisSize, omega, dim=1):
    h_HO = np.zeros((basisSize,basisSize))
    for m in range(basisSize):
        h_HO[m,m] = m + dim/2
    return h_HO * omega

def evaluate_h_DW(basisSize, omega):
    h_DW = calculate_h_plus(basisSize, omega) + \
           calculate_h_HO(basisSize, omega)
    eigval, eigvec = np.linalg.eigh(h_DW)
    return eigvec
    
def HO_SPF(n, x, omega):
    hermite_n = hermite_coefficients(n)
    summer = 0
    for i in range(len(hermite_n)):
        summer += hermite_n[i] * x ** i
    return summer * np.exp(-0.5 * omega * x * x)
    
def DW_SPF(n, L, x, omega):
    C = evaluate_h_DW(L, omega)
    summer = 0
    for lamb in range(L):
        summer += C[n,lamb] * HO_SPF(lamb, x, omega)
    return summer
    
if __name__ == '__main__':
    x = np.linspace(-5,5,1000)
    SPF = DW_SPF(3,4,x,1)
    
    import matplotlib.pyplot as plt
    plt.plot(x,abs(SPF)**2)
    plt.show()

