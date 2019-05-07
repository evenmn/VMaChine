#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system);
    double  calculateDistanceMatrixElement  (const unsigned int i, const unsigned int j);
    double  calculateRadialVectorElement    (const unsigned int particle);
    void    calculateDistanceMatrix         ();
    void    calculateRadialVector           ();
    void    setupInitialState               ();
};

