#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState
{
public:
    RandomNormal(System *system);
    double calculateDistanceMatrixElement(arma::uword i, arma::uword j);
    void calculateDistanceMatrix();
    double calculateRadialVectorElement(arma::uword particle);
    void calculateRadialVector();
    void setupInitialState();

private:
    double m_maxRadius = 1;
    double m_omega = 1;
};
