#pragma once
#include "metropolis.h"

class BruteForce : public Metropolis
{
public:
    BruteForce(System *system);
    bool acceptMove();
};
