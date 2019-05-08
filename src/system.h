#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

class System {
public:
    void runIterations                  (const unsigned int numberOfIterations);
    void printToTerminal                (const unsigned int numberOfIterations);

    void runMetropolisCycles            ();
    void checkingConvergence            ();

    void setNumberOfParticles           (const unsigned int numberOfParticles);
    void setNumberOfDimensions          (const unsigned short numberOfDimensions);
    void setNumberOfHiddenNodes         (const unsigned int numberOfHiddenNodes);
    void setNumberOfFreeDimensions      ();
    void setNumberOfMetropolisSteps     (const unsigned long steps);
    void setMaxNumberOfParametersPerElement ();
    void setNumberOfWaveFunctionElements (const unsigned long numberOfWaveFunctionElements);
    void setStepLength                  (const double stepLength);
    void setEquilibrationFraction       (const double equilibrationFraction);
    void setFrequency                   (const double omega);
    void setAtomicNumber                (const unsigned int Z);
    void setWidth                       (const double sigma);
    void setLearningRate                (const double learningRate);
    void setPath                        (const std::string path);
    void setGradients                   ();
    void setGlobalArraysToCalculate     ();

    void setInteraction                 (const bool interaction);
    void setParameterPrintingTools      (bool printParametersToFile);
    void setEnergyPrintingTools         (bool printEnergyFile, bool printInstantEnergyFile);
    void setConvergenceTools            (bool checkConvergence, unsigned int numberOfEnergies, double tolerance);
    void setDynamicStepTools            (bool applyAdaptiveSteps, unsigned int rangeOfDynamicSteps, unsigned int additionalSteps, unsigned int additionalStepsLastIteration);
    void setDensityTools                (bool computeOneBodyDensity, bool computeTwoBodyDensity, int numberOfBins, double maxRadius);
    void setMPITools                    (int myRank, int numberOfProcesses);

    unsigned long  adaptiveSteps        ();

    void setHamiltonian                 (class Hamiltonian* hamiltonian);
    void setBasis                       (class Basis* basis);
    void setInitialState                (class InitialState* initialState);
    void setInitialWeights              (class InitialWeights* initialWeights);
    void setMetropolis                  (class Metropolis* metropolis);
    void setOptimization                (class Optimization* optimization);
    void setRandomNumberGenerator       (class RandomNumberGenerator* randomNumberGenerator);
    void setWaveFunctionElements        (std::vector<class WaveFunction*> waveFunctionElements);

    class WaveFunction*                 getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*                  getHamiltonian()             { return m_hamiltonian; }
    class Basis*                        getBasis()                   { return m_basis; }
    class Sampler*                      getSampler()                 { return m_sampler; }
    class Optimization*                 getOptimization()            { return m_optimization; }
    class InitialWeights*               getInitialWeights()          { return m_initialWeights; }
    class InitialState*                 getInitialState()            { return m_initialState; }
    class RandomNumberGenerator*        getRandomNumberGenerator()   { return m_randomNumberGenerator; }

    int                                 getNumberOfProcesses()       { return m_numberOfProcesses; }
    int                                 getNumberOfBins()            { return m_numberOfBins; }
    int                                 getRank()                    { return m_myRank; }
    unsigned short                      getNumberOfDimensions()      { return m_numberOfDimensions; }
    unsigned int                        getNumberOfParticles()       { return m_numberOfParticles; }
    unsigned int                        getNumberOfHiddenNodes()     { return m_numberOfHiddenNodes; }
    unsigned int                        getNumberOfFreeDimensions()  { return m_numberOfFreeDimensions; }
    unsigned int                        getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    unsigned int                        getMaxNumberOfParametersPerElement() { return m_maxNumberOfParametersPerElement; }
    unsigned int                        getAtomicNumber()            { return m_Z; }

    unsigned long                       getNumberOfWaveFunctionElements() { return m_numberOfWaveFunctionElements; }
    unsigned long                       getTotalNumberOfStepsWOEqui(){ return m_totalNumberOfStepsWOEqui; }
    unsigned long                       getNumberOfStepsWOEqui()     { return m_numberOfStepsWOEqui; }
    unsigned long                       getTotalNumberOfStepsWEqui() { return m_totalNumberOfStepsWEqui; }
    unsigned long                       getNumberOfStepsWEqui()      { return m_numberOfStepsWEqui; }
    unsigned long                       getInitialTotalNumberOfStepsWOEqui() { return m_initialTotalNumberOfStepsWOEqui; }
    unsigned int                        getTotalNumberOfEquilibriationSteps() { return m_totalNumberOfEquilibriationSteps; }
    unsigned int                        getnumberOfEquilibriationSteps() { return m_numberOfEquilibriationSteps; }

    double                              getMaxRadius()               { return m_maxRadius; }
    double                              getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                              getFrequency()               { return m_omega; }
    double                              getWidth()                   { return m_sigma; }
    double                              getLearningRate()            { return m_eta; }
    double                              getStepLength()              { return m_stepLength; }
    bool                                getInteraction()             { return m_interaction; }
    bool                                getDensity()                 { return m_computeOneBodyDensity; }
    bool                                computeTwoBodyDensity()      { return m_computeTwoBodyDensity; }
    bool                                getPrintEnergy()             { return m_printEnergyFile; }
    bool                                getPrintParametersToFile()   { return m_printParametersToFile; }
    bool                                getPrintInstantEnergy()      { return m_printInstantEnergyFile; }
    bool                                getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool                                getCalculateRadialVector()   { return m_calculateRadialVector; }
    Eigen::VectorXd                     getPositions()               { return m_positions; }
    Eigen::VectorXd                     getRadialVector()            { return m_radialVector; }
    Eigen::MatrixXd                     getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::MatrixXd                     getWeights()                 { return m_parameters; }
    std::string                         getPath()                    { return m_path; }
    std::vector<class WaveFunction*>    getWaveFunctionElements()    { return m_waveFunctionElements; }

    void                                updateAllArrays              (const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const unsigned int changedCoord);
    void                                resetAllArrays               ();
    void                                updateAllParameters          (const Eigen::MatrixXd parameters);
    double                              evaluateWaveFunctionRatio    ();
    double                              getKineticEnergy             ();
    Eigen::MatrixXd                     getAllInstantGradients       ();
    std::string                         getAllLabels                 ();


private:
    int                                 m_numberOfProcesses         = 1;
    int                                 m_myRank                    = 0;
    int                                 m_numberOfBins              = 1;
    unsigned short                      m_numberOfDimensions        = 0;
    unsigned int                        m_numberOfHiddenNodes       = 0;
    unsigned int                        m_numberOfParticles         = 0;
    unsigned int                        m_numberOfFreeDimensions    = 0;
    unsigned int                        m_maxNumberOfParametersPerElement = 0;
    unsigned int                        m_totalNumberOfParameters   = 0;
    unsigned int                        m_Z                         = 1;
    unsigned int                        m_rangeOfDynamicSteps       = 10;
    unsigned int                        m_additionalSteps           = 4;
    unsigned int                        m_additionalStepsLastIteration = 8;
    unsigned int                        m_lastIteration             = 1;
    unsigned int                        m_numberOfEnergies          = 0;
    unsigned int                        m_iter                      = 0;
    unsigned long                       m_numberOfWaveFunctionElements = 0;

    unsigned long                       m_totalNumberOfStepsWOEqui  = 0;
    unsigned long                       m_numberOfStepsWOEqui       = 0;
    unsigned long                       m_totalNumberOfStepsWEqui   = 0;
    unsigned long                       m_numberOfStepsWEqui        = 0;
    unsigned long                       m_initialTotalNumberOfStepsWOEqui = 0;
    unsigned long                       m_initialNumberOfStepsWOEqui = 0;
    unsigned int                        m_totalNumberOfEquilibriationSteps = 0;
    unsigned int                        m_numberOfEquilibriationSteps = 0;

    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    double                              m_omega                     = 1.0;
    double                              m_sigma                     = 1.0;
    double                              m_eta                       = 0.1;
    double                              m_tolerance                 = 1e-7;
    double                              m_maxRadius                 = 10;
    double                              m_totalTime                 = 0;

    bool                                m_interaction               = true;
    bool                                m_checkConvergence          = true;
    bool                                m_applyAdaptiveSteps        = true;
    bool                                m_computeOneBodyDensity     = true;
    bool                                m_computeTwoBodyDensity     = true;
    bool                                m_printEnergyFile           = true;
    bool                                m_printInstantEnergyFile    = true;
    bool                                m_printParametersToFile     = true;
    bool                                m_calculateDistanceMatrix   = true;
    bool                                m_calculateRadialVector     = true;

    std::string                         m_path                      = "data/";

    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class Basis*                        m_basis                     = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    class Metropolis*                   m_metropolis                = nullptr;
    class Optimization*                 m_optimization              = nullptr;
    class RandomNumberGenerator*        m_randomNumberGenerator     = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionElements;

    Eigen::VectorXd                     m_positions;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
    Eigen::VectorXd                     m_energies;
};

