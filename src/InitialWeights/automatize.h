#pragma once
#include <memory>

#include "initialweights.h"

class Automatize : public InitialWeights
{
public:
    Automatize(System *system);
    std::string getLabel() { return m_method->getLabel(); }
    void setupInitialWeights();

    std::string generateFileName(std::string name, std::string extension);
    Eigen::MatrixXd getParameters();

private:
    double m_factor = 0;
    std::string m_trialWaveFunction;
    std::string m_hamiltonian;

    int m_initialTotalStepsWOEqui = 0;
    bool m_interaction = 0;
    std::string m_path = "No path";

    std::unique_ptr<class InitialWeights> m_method = nullptr;
};
