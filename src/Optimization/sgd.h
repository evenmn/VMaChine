#pragma once
#include "optimization.h"

class SGD : public Optimization
{
public:
    SGD(System *system, const double gamma = 0, const double monotonicExp = 0);
    int getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel() { return m_label; }

    void initialize();
    Eigen::MatrixXd updateParameters();

private:
    int m_step = 0;
    int m_numberOfBatches = 1;
    double m_gamma = 0;
    double m_monotonicExp = 0;
    Eigen::MatrixXd m_v;
    std::string m_label = "stochastic gradient descent";
};
