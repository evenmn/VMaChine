#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState
{
public:
    RandomNormal(System *system);
    void setupInitialState();

private:
    double m_maxRadius = 1;
    double m_omega = 1;
};
