#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState {
public:
    RandomNormal(System* system);
    double  calculateDistanceMatrixElement  (const unsigned int i, const unsigned int j);
    double  calculateRadialVectorElement    (const unsigned int particle);
    void    calculateDistanceMatrix         ();
    void    calculateRadialVector           ();
    void    setupInitialState               ();
};
