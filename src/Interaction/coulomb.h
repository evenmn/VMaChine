#pragma once
#include "interaction.h"
#include "../system.h"

class Coulomb : public Interaction
{
public:
    Coulomb(System *system);
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double getInteractionEnergy();

private:
    int m_globalArrayNeed = 1;
    std::string m_label = "Coulomb";

    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    bool m_screening = false;

    double m_screeningStrength = 100;
    double m_dsl = 100;

    Eigen::MatrixXd m_distanceMatrix;

    Eigen::Map<Eigen::VectorXd> flat(Eigen::MatrixXd A);

    // class System *m_system = nullptr;
};
