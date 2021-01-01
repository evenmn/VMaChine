#include "nointeraction.h"
#include <cassert>

NoInteraction::NoInteraction(System *system)
    : Interaction(system)
{}

void NoInteraction::initialize() {}

double NoInteraction::getInteractionEnergy()
{
    return 0;
}
