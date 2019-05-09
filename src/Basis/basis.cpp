#include "basis.h"

Basis::Basis(System *system) {
}

long long Basis::factorial(const int n) {
    return (n <= 1) ? 1 : factorial(n - 1) * n;
}

int Basis::factorialDifference(const int high, const int low) {
    return (high <= 1 || high <= low) ? 1 : factorialDifference(high - 1, low) * high;
}

int Basis::binomial(const int n, const int p) {
    //Binomial coefficients, equal to magic numbers
    //return int(factorial(n+p)/(factorial(n)*factorial(p)));
    return factorialDifference(n+p, n) / factorial(p);
}

Basis::~Basis() {};
