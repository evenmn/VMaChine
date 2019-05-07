#include "basis.h"

Basis::Basis(System *system) {
}

unsigned long long Basis::factorial(const unsigned int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

unsigned int Basis::binomial(const unsigned int n, const unsigned int p) {
    //Binomial coefficients, equal to magic numbers
    return unsigned(factorial(n+p)/(factorial(n)*factorial(p)));
}

Basis::~Basis() {};
