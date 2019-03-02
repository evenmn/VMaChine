#include "hermite.h"
#include "../system.h"
#include <iostream>

Hermite::Hermite()  :
    Basis() {
}

double Hermite::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    //std::cout << n << std::endl;
    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 2*x;
    }
    else {
        return 2*x*evaluate(x,n-1)-2*(n-1)*evaluate(x,n-2);
    }
}

double Hermite::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    if(n == 0) {
        return 0;
    }
    else {
        return 2*n*evaluate(x,n-1);
    }
}
