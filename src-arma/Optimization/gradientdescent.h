#pragma once
#include "optimization.h"

class GradientDescent : public Optimization
{
public:
    GradientDescent(System *system, const double gamma = 0, const double monotonicExp = 0);
    arma::uword getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel() { return m_label; }

    void initialize();
    arma::mat updateParameters();

private:
    arma::uword m_step = 0;
    arma::uword m_numberOfBatches = 1;
    double m_gamma = 0;
    double m_monotonicExp = 0;
    arma::mat m_v;

    std::string m_label = "GD";
};
