#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState
{
public:
    RandomUniform(System *system);
    double calculateDistanceMatrixElement(int i, int j);
    void calculateDistanceMatrix();
    double calculateRadialVectorElement(int particle);
    void calculateRadialVector();
    void setupInitialState();
};
