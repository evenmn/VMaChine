#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator {
public:
    MersenneTwister();
    unsigned int nextInt     (unsigned int upperLimit);
    double       nextDouble  ();
    double       nextGaussian(double mean, double standardDeviation);
};
