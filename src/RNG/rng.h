#pragma once
#include <Eigen/Dense>

class RandomNumberGenerator {
public:
    RandomNumberGenerator();
    virtual unsigned int    nextInt     (unsigned int upperLimit) = 0;
    virtual double          nextDouble  () = 0;
    virtual double          nextGaussian(double mean, double standardDeviation) = 0;

    virtual ~RandomNumberGenerator() = 0;
};
