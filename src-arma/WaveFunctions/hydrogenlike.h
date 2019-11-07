#pragma once
#include "wavefunction.h"

class HydrogenLike : public WaveFunction
{
public:
    HydrogenLike(class System *system);
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
    arma::uword m_globalArrayNeed = 2;
    arma::uword m_elementNumber = 0;
    arma::uword m_Z = 1;

    double m_alpha = 0;

    arma::vec m_positions;
    arma::vec m_positionsOld;
    arma::vec m_radialVector;
    arma::vec m_radialVectorOld;
    arma::vec m_gradients;

    double m_probabilityRatio;
    double m_probabilityRatioOld;

    std::string m_label = "hydrogenlike";
};
