#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

class System {
public:
    void runMetropolisSteps             (const int numberOfIterations);
    void setNumberOfParticles           (const int numberOfParticles);
    void setNumberOfDimensions          (const int numberOfDimensions);
    void setNumberOfHiddenNodes         (const int numberOfHiddenNodes);
    void setNumberOfOrbitals            ();
    void setNumberOfFreeDimensions      ();
    void setTotalNumberOfSteps          ();
    void setNumberOfMetropolisSteps     (const int steps);
    void setMaxNumberOfParametersPerElement (const int maxNumberOfParametersPerElement);
    void setNumberOfWaveFunctionElements (const int numberOfWaveFunctionElements);
    void setStepLength                  (const double stepLength);
    void setEquilibrationFraction       (const double equilibrationFraction);
    void setFrequency                   (const double omega);
    void setWidth                       (const double sigma);
    void setLearningRate                (const double eta);
    void setInteraction                 (const bool interaction);
    void setGradients                   ();
    void setHamiltonian                 (class Hamiltonian* hamiltonian);
    void setInitialState                (class InitialState* initialState);
    void setInitialWeights              (class InitialWeights* initialWeights);
    void setMetropolis                  (class Metropolis* metropolis);
    void setOptimization                (class Optimization* optimization);
    void setRandomNumberGenerator       (class RandomNumberGenerator* randomnumbergenerator);
    void setWaveFunction                (std::vector<class WaveFunction*> waveFunctionVector);

    class WaveFunction*                 getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*                  getHamiltonian()             { return m_hamiltonian; }
    class Sampler*                      getSampler()                 { return m_sampler; }
    class Optimization*                 getOptimization()            { return m_optimization; }
    class InitialWeights*               getInitialWeights()          { return m_initialWeights; }
    class RandomNumberGenerator*        getRandomNumberGenerator()   { return m_randomnumbergenerator; }
    int                                 getNumberOfParticles()       { return m_numberOfParticles; }
    int                                 getNumberOfDimensions()      { return m_numberOfDimensions; }
    int                                 getNumberOfHiddenNodes()     { return m_numberOfHiddenNodes; }
    int                                 getNumberOfOrbitals()        { return m_numberOfOrbitals; }
    int                                 getNumberOfFreeDimensions()  { return m_numberOfFreeDimensions; }
    int                                 getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    int                                 getTotalNumberOfSteps()      { return m_totalNumberOfSteps; }
    int                                 getMaxNumberOfParametersPerElement() { return m_maxNumberOfParametersPerElement; }
    int                                 getNumberOfWaveFunctionElements() { return m_numberOfWaveFunctionElements; }
    double                              getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                              getFrequency()               { return m_omega; }
    double                              getWidth()                   { return m_sigma; }
    double                              getLearningRate()            { return m_eta; }
    double                              getStepLength()              { return m_stepLength; }
    bool                                getInteraction()             { return m_interaction; }
    Eigen::VectorXd                     getParticles()               { return m_positions; }
    Eigen::VectorXd                     getPositions()               { return m_positions; }
    Eigen::MatrixXd                     getWeights()                 { return m_parameters; }
    Eigen::MatrixXd                     getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::VectorXd                     getRadialVector()            { return m_radialVector; }
    std::vector<class WaveFunction*>    getWaveFunctionElements()    { return m_waveFunctionVector; }

    void            updateAllArrays          (const Eigen::VectorXd positions, const int pRand);
    void            resetAllArrays           ();
    void            updateAllParameters      (const Eigen::MatrixXd parameters);
    double          evaluateWaveFunction     ();
    double          evaluateWaveFunctionSqrd ();
    double          getKineticEnergy         ();
    std::string     generate_filename        (std::string name, const std::string extension);

private:
    int                                 m_numberOfHiddenNodes       = 0;
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfFreeDimensions    = 0;
    int                                 m_numberOfMetropolisSteps   = 0;
    int                                 m_numberOfOrbitals          = 0;
    int                                 m_numberOfWaveFunctionElements = 0;
    int                                 m_maxNumberOfParametersPerElement = 0;
    int                                 m_totalNumberOfSteps        = 0;
    bool                                m_interaction               = false;
    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    double                              m_omega                     = 1.0;
    double                              m_sigma                     = 1.0;
    double                              m_eta                       = 0.1;
    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    class Metropolis*                   m_metropolis                = nullptr;
    class Optimization*                 m_optimization              = nullptr;
    class RandomNumberGenerator*        m_randomnumbergenerator     = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
    Eigen::VectorXd                     m_positions;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
};

