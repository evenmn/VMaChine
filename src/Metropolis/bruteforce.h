#pragma once
#include "metropolis.h"
#include <Eigen/Dense>

class BruteForce : public Metropolis {
public:
    BruteForce(System* system);
    bool acceptMove();

    double calculateDistanceMatrixElement(int i, int j);
    void   calculateDistanceMatrixCross(int particle);
    double calculateRadialVectorElement(int particle);
};
