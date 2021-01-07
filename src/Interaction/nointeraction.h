#pragma once
#include "interaction.h"

class NoInteraction : public Interaction
{
public:
    NoInteraction(System *system);
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double getInteractionEnergy();

private:
    int m_globalArrayNeed = 0;
    std::string m_label = "no interaction";
};
