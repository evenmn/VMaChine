#pragma once
#include "optimization.h"

class ASGD : public Optimization
{
public:
    ASGD(System *system, const double gamma);
    arma::uword getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel() { return m_label; }

    void initialize();
    arma::mat updateParameters();

private:
    double m_A = 20;
    double m_gamma = 0;
    double m_t = m_A;
    arma::uword m_step = 0;
    arma::uword m_numberOfBatches = 10;
    arma::mat m_v;

    std::string m_label = "ASGD";
};
