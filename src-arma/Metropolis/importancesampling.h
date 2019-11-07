#pragma once
#include "metropolis.h"
#include <vector>

class ImportanceSampling : public Metropolis
{
public:
    ImportanceSampling(System *system);
    void initialize();
    bool acceptMove();

    void initializeQuantumForce();
    double QuantumForce(const arma::uword i);
    double GreenRatio(const arma::uword i);

private:
    arma::vec m_quantumForceOld;
    arma::vec m_quantumForceNew;
    std::vector<class WaveFunction *> m_waveFunctionVector;

    double m_dtD = 1;
    double m_dx = 0;
    double m_sqrtStep = 0;
};
