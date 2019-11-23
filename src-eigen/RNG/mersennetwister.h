#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int nextInt(int upperLimit);
    double nextDouble();
    double nextGaussian(double mean, double standardDeviation);
};
