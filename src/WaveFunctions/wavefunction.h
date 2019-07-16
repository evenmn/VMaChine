#pragma once
#include "../Eigen/Dense"
#include <iostream>
#include <string>

class WaveFunction {

public:
    WaveFunction(class System *system);
    virtual int             getNumberOfParameters       () = 0;
    virtual int             getGlobalArrayNeed          () = 0;
    virtual std::string     getLabel                    () = 0;

    virtual void            updateParameters            (const Eigen::MatrixXd parameters) = 0;
    virtual void            initializeArrays            (const Eigen::VectorXd positions, \
                                                         const Eigen::VectorXd radialVector, \
                                                         const Eigen::MatrixXd distanceMatrix) = 0;
    virtual void            updateArrays                (const Eigen::VectorXd positions, \
                                                         const Eigen::VectorXd radialVector, \
                                                         const Eigen::MatrixXd distanceMatrix, \
                                                         const int changedCoord) = 0;
    virtual void            setConstants                (const int elementNumber) = 0;
    virtual void            setArrays                   () = 0;
    virtual void            resetArrays                 () = 0;
    virtual double          evaluateRatio               () = 0;
    virtual double          computeGradient             (const int k) = 0;
    virtual double          computeLaplacian            () = 0;
    virtual Eigen::VectorXd computeParameterGradient    () = 0;

    virtual ~WaveFunction() = 0;

    Eigen::Map<Eigen::VectorXd> flatten(Eigen::MatrixXd A);
    Eigen::Map<Eigen::MatrixXd> reshape(Eigen::VectorXd A, const Eigen::Index m, const Eigen::Index n);
    Eigen::Map<Eigen::MatrixXd> square(Eigen::VectorXd A);

protected:
    int     m_numberOfParticles         = 0;
    int     m_numberOfDimensions        = 0;
    int     m_numberOfFreeDimensions    = 0;
    int     m_maxParameters             = 0;
    class System* m_system = nullptr;
};

