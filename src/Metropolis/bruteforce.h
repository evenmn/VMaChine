#pragma once
#include "metropolis.h"
#include <Eigen/Dense>

class BruteForce : public Metropolis {
public:
    BruteForce(System* system);
    bool acceptMove();
};
