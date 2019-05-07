#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <string>

class WaveFunction {

public:
    WaveFunction(class System *system);
    virtual unsigned int    getNumberOfParameters       ()                                   = 0;
    virtual unsigned short  getGlobalArrayNeed          ()                                   = 0;
    virtual std::string     getLabel                    ()                                   = 0;

    virtual void            updateArrays                (const Eigen::VectorXd positions,      \
                                                         const Eigen::VectorXd radialVector,   \
                                                         const Eigen::MatrixXd distanceMatrix, \
                                                         const unsigned int changedCoord)    = 0;
    virtual void            setArrays                   ()                                   = 0;
    virtual void            resetArrays                 ()                                   = 0;
    virtual void            initializeArrays            (const Eigen::VectorXd positions,      \
                                                         const Eigen::VectorXd radialVector,   \
                                                         const Eigen::MatrixXd distanceMatrix) = 0;
    virtual void            updateParameters            (const Eigen::MatrixXd parameters,     \
                                                         const unsigned short elementNumber) = 0;
    virtual double          evaluateRatio               ()                                   = 0;
    virtual double          computeGradient             (const unsigned int k)               = 0;
    virtual double          computeLaplacian            ()                                   = 0;
    virtual Eigen::VectorXd computeParameterGradient    ()                                   = 0;

    virtual ~WaveFunction() = 0;

protected:
    unsigned short          m_numberOfDimensions                = 0;
    unsigned int            m_numberOfParticles                 = 0;
    unsigned int            m_numberOfFreeDimensions            = 0;
    unsigned int            m_maxNumberOfParametersPerElement   = 0;

    class System*           m_system                            = nullptr;
};

