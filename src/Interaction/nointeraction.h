#pragma once
#include "interaction.h"

class NoInteraction : public Interaction
{
public:
    NoInteraction(System *system);
    int getGlobalArrayNeed() override { return m_globalArrayNeed; }
    std::string getLabel() override { return m_label; }

    void initialize() override;
    double getInteractionEnergy() override;

private:
    int m_globalArrayNeed = 0;
    std::string m_label = "no interaction";
};
