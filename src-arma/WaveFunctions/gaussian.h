#pragma once
#include "wavefunction.h"

class Gaussian : public WaveFunction
{
public:
    Gaussian(class System *system);
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

    void updateProbabilityRatio(arma::uword changedCoord);

private:
    arma::uword m_numberOfParameters = 1;
    arma::uword m_globalArrayNeed = 0;
    arma::uword m_elementNumber = 0;

    double m_omega = 0;
    double m_alpha = 0;
    double m_omegalpha = 0;
    double m_probabilityRatio = 0;
    double m_probabilityRatioOld = 0;
    arma::vec m_positions;
    arma::vec m_positionsOld;

    arma::vec m_gradients;

    std::string m_label = "gaussian";
};
