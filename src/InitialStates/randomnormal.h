#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState {
public:
    RandomNormal(System* system);
    double  calculateDistanceMatrixElement  (int i, int j);
    void    calculateDistanceMatrix         ();
    double  calculateRadialVectorElement    (int particle);
    void    calculateRadialVector           ();
    void    setupInitialState               ();
};
