#pragma once
#include "metropolis.h"
#include <Eigen/Dense>

class BruteForce : public Metropolis {
public:
    BruteForce(System* system);
    bool acceptMove();

    double calculateDistanceMatrixElement(const unsigned int i, const unsigned int j);
    void   calculateDistanceMatrixCross(const unsigned int particle);
    double calculateRadialVectorElement(const unsigned int particle);
};
