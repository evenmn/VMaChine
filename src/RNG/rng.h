#pragma once

class RandomNumberGenerator {
public:
    RandomNumberGenerator();
    virtual int    nextInt(int upperLimit) = 0;
    virtual double nextDouble() = 0;
    virtual double nextGaussian(double mean, double standardDeviation) = 0;

    virtual ~RandomNumberGenerator() = 0;
};
