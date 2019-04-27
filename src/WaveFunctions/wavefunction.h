#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <string>

class WaveFunction {

public:
    WaveFunction(class System *system);
    virtual int             getNumberOfParameters           ()                                      = 0;
    virtual int             getGlobalArrayNeed              ()                                      = 0;
    virtual std::string     getLabel                        ()                                      = 0;

    virtual void            updateArrays                    (Eigen::VectorXd positions, \
                                                             Eigen::VectorXd radialVector, \
                                                             Eigen::MatrixXd distanceMatrix, int pRand)  = 0;
    virtual void            setArrays                       ()                                      = 0;
    virtual void            resetArrays                     ()                                      = 0;
    virtual void            initializeArrays                (Eigen::VectorXd positions, \
                                                             Eigen::VectorXd radialVector, \
                                                             Eigen::MatrixXd distanceMatrix)        = 0;
    virtual void            updateParameters                (Eigen::MatrixXd parameters, const int elementNumber)            = 0;
    virtual double          evaluateRatio                   ()                                      = 0;
    virtual double          computeGradient                 (int k)                                 = 0;
    virtual double          computeLaplacian                ()                                      = 0;
    virtual Eigen::VectorXd computeParameterGradient        ()                                      = 0;

    virtual ~WaveFunction() = 0;

protected:
    int     m_numberOfParticles                 = 0;
    int     m_numberOfDimensions                = 0;
    int     m_numberOfFreeDimensions            = 0;
    int     m_maxNumberOfParametersPerElement   = 0;
    class System* m_system = nullptr;
};

