#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator {
public:
    MersenneTwister(class System* system);
    int    nextInt(int upperLimit);
    double nextDouble();
    double nextGaussian(double mean, double standardDeviation);
};
