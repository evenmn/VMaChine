#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

class System {
public:
    void runIterations                  (const int numberOfIterations);
    void runMetropolisCycles            (int numberOfSteps, int totalNumberOfSteps, int iter);
    void printToTerminal                (int numberOfIterations, int iter, double time);
    void checkingConvergence            (int iter);
    int  dynamicSteps                   (int iter);
    void setNumberOfParticles           (const int numberOfParticles);
    void setNumberOfDimensions          (const int numberOfDimensions);
    void setNumberOfHiddenNodes         (const int numberOfHiddenNodes);
    void setNumberOfFreeDimensions      ();
    void setTotalNumberOfSteps          ();
    void setNumberOfMetropolisSteps     (const int steps);
    void setMaxNumberOfParametersPerElement ();
    void setNumberOfWaveFunctionElements (const int numberOfWaveFunctionElements);
    void setStepLength                  (const double stepLength);
    void setEquilibrationFraction       (const double equilibrationFraction);
    void setFrequency                   (const double omega);
    void setAtomicNumber                (const int Z);
    void setWidth                       (const double sigma);
    void setLearningRate                (const double eta);
    void setGradients                   ();
    void setGlobalArraysToCalculate     ();

    void setInteraction                 (const bool interaction);
    void setConvergenceTools            (bool checkConvergence, int numberOfEnergies, double tolerance);
    void setDynamicStepTools            (bool applyDynamicSteps, int rangeOfDynamicSteps, int additionalSteps, int additionalStepsLastIteration);
    void setDensityTools                (bool computeDensity, int numberOfBins, double maxRadius);
    void setEnergyPrintingTools         (bool printEnergyFile, bool printInstantEnergyFile);

    void setHamiltonian                 (class Hamiltonian* hamiltonian);
    void setBasis                       (class Basis* basis);
    void setInitialState                (class InitialState* initialState);
    void setInitialWeights              (class InitialWeights* initialWeights);
    void setMetropolis                  (class Metropolis* metropolis);
    void setOptimization                (class Optimization* optimization);
    void setRandomNumberGenerator       (class RandomNumberGenerator* randomnumbergenerator);
    void setWaveFunctionElements        (std::vector<class WaveFunction*> waveFunctionElements);

    class WaveFunction*                 getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*                  getHamiltonian()             { return m_hamiltonian; }
    class Basis*                        getBasis()                   { return m_basis; }
    class Sampler*                      getSampler()                 { return m_sampler; }
    class Optimization*                 getOptimization()            { return m_optimization; }
    class InitialWeights*               getInitialWeights()          { return m_initialWeights; }
    class InitialState*                 getInitialState()            { return m_initialState; }
    class RandomNumberGenerator*        getRandomNumberGenerator()   { return m_randomnumbergenerator; }
    int                                 getNumberOfParticles()       { return m_numberOfParticles; }
    int                                 getNumberOfDimensions()      { return m_numberOfDimensions; }
    int                                 getNumberOfHiddenNodes()     { return m_numberOfHiddenNodes; }
    int                                 getNumberOfFreeDimensions()  { return m_numberOfFreeDimensions; }
    int                                 getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    int                                 getTotalNumberOfSteps()      { return m_totalNumberOfSteps; }
    int                                 getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    int                                 getMaxNumberOfParametersPerElement() { return m_maxNumberOfParametersPerElement; }
    int                                 getNumberOfWaveFunctionElements() { return m_numberOfWaveFunctionElements; }
    int                                 getAtomicNumber()            { return m_Z; }
    int                                 getNumberOfBins()            { return m_numberOfBins; }
    double                              getMaxRadius()               { return m_maxRadius; }
    double                              getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                              getFrequency()               { return m_omega; }
    double                              getWidth()                   { return m_sigma; }
    double                              getLearningRate()            { return m_eta; }
    double                              getStepLength()              { return m_stepLength; }
    bool                                getInteraction()             { return m_interaction; }
    bool                                getDensity()                 { return m_computeDensity; }
    bool                                getPrintEnergy()             { return m_printEnergyFile; }
    bool                                getPrintInstantEnergy()      { return m_printInstantEnergyFile; }
    bool                                getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool                                getCalculateRadialVector()   { return m_calculateRadialVector; }
    Eigen::VectorXd                     getPositions()               { return m_positions; }
    Eigen::VectorXd                     getRadialVector()            { return m_radialVector; }
    Eigen::MatrixXd                     getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::MatrixXd                     getWeights()                 { return m_parameters; }
    std::vector<class WaveFunction*>    getWaveFunctionElements()    { return m_waveFunctionElements; }

    void                                updateAllArrays              (const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord);
    void                                resetAllArrays               ();
    void                                updateAllParameters          (const Eigen::MatrixXd parameters);
    double                              evaluateWaveFunctionRatio    ();
    double                              getKineticEnergy             ();
    Eigen::MatrixXd                     getAllInstantGradients       ();
    std::string                         getAllLabels                 ();


private:
    int                                 m_numberOfHiddenNodes       = 0;
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfFreeDimensions    = 0;
    int                                 m_numberOfMetropolisSteps   = 0;
    int                                 m_numberOfWaveFunctionElements = 0;
    int                                 m_maxNumberOfParametersPerElement = 0;
    int                                 m_totalNumberOfSteps        = 0;
    int                                 m_totalNumberOfParameters   = 0;
    int                                 m_Z                         = 1;
    int                                 m_rangeOfDynamicSteps       = 10;
    int                                 m_additionalSteps           = 4;
    int                                 m_additionalStepsLastIteration = 8;
    int                                 m_lastIteration             = 1;
    int                                 m_numberOfEnergies          = 0;
    int                                 m_numberOfBins              = 1;

    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    double                              m_omega                     = 1.0;
    double                              m_sigma                     = 1.0;
    double                              m_eta                       = 0.1;
    double                              m_tolerance                 = 1e-7;
    double                              m_maxRadius                 = 10;

    bool                                m_interaction               = false;
    bool                                m_checkConvergence          = false;
    bool                                m_applyDynamicSteps         = false;
    bool                                m_computeDensity            = true;
    bool                                m_printEnergyFile           = true;
    bool                                m_printInstantEnergyFile    = true;
    bool                                m_calculateDistanceMatrix   = false;
    bool                                m_calculateRadialVector     = false;

    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class Basis*                        m_basis                     = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    class Metropolis*                   m_metropolis                = nullptr;
    class Optimization*                 m_optimization              = nullptr;
    class RandomNumberGenerator*        m_randomnumbergenerator     = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionElements;

    Eigen::VectorXd                     m_positions;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
    Eigen::VectorXd                     m_energies;
};

