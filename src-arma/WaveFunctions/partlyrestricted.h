#pragma once
#include "wavefunction.h"

class PartlyRestricted : public WaveFunction
{
public:
    PartlyRestricted(class System *system);
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

    void calculateProbabilityRatio(arma::uword changedCoord);

private:
    arma::uword m_numberOfParameters = 1;
    arma::uword m_globalArrayNeed = 0;
    arma::uword m_elementNumber = 0;
    double m_xCxOld;
    double m_xCx;
    double m_probabilityRatio;
    double m_probabilityRatioOld;

    arma::vec m_positions;
    arma::vec m_positionsOld;
    arma::mat m_c;
    arma::vec m_gradients;

    std::string m_label = "partlyrestricted";
};
