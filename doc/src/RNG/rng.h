#pragma once
#include <Eigen/Dense>

class RandomNumberGenerator {
public:
    RandomNumberGenerator(class System* system);
    virtual int    nextInt(int upperLimit) = 0;
    virtual double nextDouble() = 0;
    virtual double nextGaussian(double mean, double standardDeviation) = 0;

    virtual ~RandomNumberGenerator() = 0;

protected:
    class System* m_system = nullptr;
};
