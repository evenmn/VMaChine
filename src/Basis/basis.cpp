#include "basis.h"

Basis::Basis(System *system) {
}

int Basis::factorial(const int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double Basis::binomial(const int n, const int p) {
    //Binomial coefficients, equal to magic numbers
    return factorial(n+p)/(factorial(n)*factorial(p));
}

Basis::~Basis() {};
