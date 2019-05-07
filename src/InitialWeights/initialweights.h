#pragma once
#include <Eigen/Dense>

class InitialWeights {
public:
    InitialWeights(class System* system);
    virtual void    setupInitialWeights() = 0;
    Eigen::MatrixXd getWeights() { return m_parameters; }

    virtual         ~InitialWeights() = 0;

protected:
    unsigned short  m_numberOfDimensions    = 0;
    unsigned int    m_numberOfParticles     = 0;
    unsigned short  m_numberOfElements      = 0;
    unsigned int    m_maxNumberOfParametersPerElement = 0;

    Eigen::MatrixXd m_parameters;

    class System* m_system = nullptr;
};
