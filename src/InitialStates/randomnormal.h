#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState
{
public:
    RandomNormal(System *system);
    std::string getLabel() { return m_label; }
    void setupInitialState();

private:
    double m_maxRadius = 1;
    double m_omega = 1;
    std::string m_label = "random normal";
};
