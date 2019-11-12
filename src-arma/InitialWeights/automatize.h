#pragma once
#include "initialweights.h"

class Automatize : public InitialWeights
{
public:
    Automatize(System *system);
    void setupInitialWeights();

    std::string generateFileName(std::string name, std::string extension);
    arma::mat getParameters();

private:
    double m_factor = 0;
    std::string m_trialWaveFunction;
    std::string m_hamiltonian;

    arma::uword m_initialTotalStepsWOEqui = 0;
    bool interaction = 0;
    std::string m_path = "No path";

    class InitialWeights *m_method = nullptr;
};