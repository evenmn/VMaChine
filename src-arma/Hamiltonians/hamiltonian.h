#pragma once
#include <armadillo>

class Hamiltonian
{
public:
    Hamiltonian(class System *system);
    virtual arma::uword getGlobalArrayNeed() = 0;
    virtual std::string getLabel() = 0;

    virtual void initialize() = 0;
    virtual double getExternalEnergy() = 0;
    virtual double computeLocalEnergy() = 0;
    virtual ~Hamiltonian() = 0;

    double getInteractionEnergy();

protected:
    arma::uword m_numberOfParticles = 0;
    arma::uword m_numberOfDimensions = 0;
    bool m_interaction = false;
    bool m_screening = false;

    double m_screeningStrength = 100;
    double m_dsl = 100;

    arma::vec m_positions;
    arma::vec m_radialVector;
    arma::mat m_distanceMatrix;
    arma::mat m_parameters;

    class System *m_system = nullptr;
};
