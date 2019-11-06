#pragma once
#include "wavefunction.h"

class RBMGaussian : public WaveFunction
{
public:
    RBMGaussian(class System *system);
    arma::uword getNumberOfParameters() { return m_numberOfParameters; }
    arma::uword getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void updateParameters(const arma::mat parameters);
    void initializeArrays(const arma::vec positions,
                          const arma::vec radialVector,
                          const arma::mat distanceMatrix);
    void updateArrays(const arma::vec positions,
                      const arma::vec radialVector,
                      const arma::mat distanceMatrix,
                      const arma::uword changedCoord);
    void setConstants(const arma::uword elementNumber);
    void setArrays();
    void resetArrays();
    double evaluateRatio();
    double computeGradient(const arma::uword k);
    double computeLaplacian();
    arma::vec computeParameterGradient();

private:
    arma::uword m_numberOfParameters = 1;
    arma::uword m_globalArrayNeed = 0;
    arma::uword m_elementNumber = 0;
    double m_omega = 1;
    double m_sigmaSqrd = 1;

    double m_probabilityRatio = 0;
    double m_probabilityRatioOld = 0;

    arma::vec m_positions;
    arma::vec m_positionsOld;
    arma::vec m_Xa;
    arma::vec m_XaOld;
    arma::vec m_a;
    arma::vec m_gradients;

    std::string m_label = "rbmgaussian";
};
