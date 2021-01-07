#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState
{
public:
    RandomUniform(System *system);
    std::string getLabel() { return m_label; }
    void setupInitialState();

private:
    std::string m_label = "random uniform";
};
