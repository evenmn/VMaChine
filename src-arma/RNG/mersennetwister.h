#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    unsigned long long nextInt(unsigned long long upperLimit);
    double nextDouble();
    double nextGaussian(double mean, double standardDeviation);
};
