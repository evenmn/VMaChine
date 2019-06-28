from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing, simplify, acos, sin, cos, pretty_print

x, y, z, r, k, theta, phi = symbols('x y z r k theta phi')

r = sqrt(x*x+y*y+z*z)
theta = acos(z/r)
R = Symbol('r')
        
def associatedLaguerre(X, p, q):
    if q==0:
        return 1
    elif q==1:
        return 1 + p - X
    elif q==2:
        return 0.5*(X*X - 2*(p+2)*X + (p+1)*(p+2))
    elif q==3:
        return (1/6.)*(-X*X*X + 3*(p+3)*X*X - 3*(p+2)*(p+3)*X + (p+1)*(p+2)*(p+3))
    else:
        return ((2*(q-1)+1+p-X)*associatedLaguerre(X,p,q-1) - (q-1+p)*associatedLaguerre(X,p,q-2))/q

def associatedLegendre(l, m):
    # Associated Legendre polynomials
    m = abs(m)
    if l==0:
        return 1
    elif l==1:
        if m==0:
            return z/r
        elif m==1:
            return -sin(theta)
    elif l==2:
        if m==0:
            return 0.5*(3*cos(theta)**2-1)
        elif m==1:
            return -3*sin(theta)*cos(theta)
        elif m==2:
            return 3*sin(theta)**2
    elif l==3:
        if m==0:
            return 0.5*cos(theta)*(5*cos(theta**2-3))
        elif m==1:
            return -1.5*sin(theta)*(5*cos(theta)**2-1)
        elif m==2:
            return 15*cos(theta)*sin(theta)**2
        elif m==3:
            return -15*sin(theta)**3
            
def cosine(m):
    # Returns cos(m*phi) with the assumption cos(phi)=z
    cosphi = x/(r*sin(theta))
    if m==0:
        return 1
    elif m==1:
        return cosphi
    else:
        return 2*cosphi*cosine(m-1)-cosine(m-2)
        
def sine(m):
    # Returns sin(m*phi) with the assumption cos(phi)=z and sin(phi)=x^2+y^2
    cosphi = x/(r*sin(theta))
    sinphi = y/(r*sin(theta))
    if m==0:
        return 0
    elif m==1:
        return sinphi
    else:
        return 2*cosphi*sine(m-1)-sine(m-2)
            
def solidHarmonics(l, m):
    if m>=0:
        return r**l*associatedLegendre(l, m) * cosine(m)
    else:
        return r**l*associatedLegendre(l, m) * sine(abs(m))
        
for n in range(1,6):
    N = n
    if n==5:
        N = 2
    for l in range(N):
        for m in range(-l,l+1):
            Phi = exp(-k*r/n)*associatedLaguerre(2*k*r/n, 2*l+1, n-l-1)*solidHarmonics(l, m)
            
            Derx = diff(Phi, x)
            Dery = diff(Phi, y)
            Derz = diff(Phi, z)
            
            Derxx = diff(Derx, x)
            Deryy = diff(Dery, y)
            Derzz = diff(Derz, z)
            
            print("n: ", n, " | l: ", l, " | m: ", m)
            print("==============================")
            print(" ")
            print("φ:")
            print(Phi.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_x φ:")
            print(Derx.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_y φ:")
            print(Dery.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_z φ:")
            print(Derz.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_xx φ:")
            print(Derxx.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_yy φ:")
            print(Deryy.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print(" ")
            print("∇_zz φ:")
            print(Derzz.factor().factor().simplify().simplify().subs(r, R).subs(r, R))
            print("\n\n\n")
