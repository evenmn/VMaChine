#pragma once
#include "metropolis.h"

class BruteForce : public Metropolis
{
public:
    BruteForce(System *system);
    void initialize();
    bool acceptMove();
    std::string getLabel() { return m_label; }

private:
    std::string m_label = "brute force Metropolis";
};
