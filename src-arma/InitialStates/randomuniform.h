#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState
{
public:
    RandomUniform(System *system);
    double calculateDistanceMatrixElement(arma::uword i, arma::uword j);
    void calculateDistanceMatrix();
    double calculateRadialVectorElement(arma::uword particle);
    void calculateRadialVector();
    void setupInitialState();
};
