#pragma once
#include <cassert>
#include <ctime>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mpi.h>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "main.h"

class System
{
    /* The system class is the primary class of the software.
     * This is where the parameter loop and Monte Carlo loop
     * are found. All communication goes from this class.
     * Additionally, this class stores all the required
     * information about the system. This includes particle
     * configurations, parameter values etc.. */
public:
    /* Constructor */
    System();

    /* Destructor */
    ~System();

    /* Call this function from main.cpp to specify the number
     * of particles used in the simulation. */
    void setNumberOfParticles(const int numberOfParticles);

    /* Call this function from main.cpp to specify the number
     * of dimensions of the system. */
    void setNumberOfDimensions(const int numberOfDimensions);

    /* Call this function from main.cpp to specify the number
     * of hidden nodes contained in the hidden layer of the
     * restricted Boltzmann machine. The default setting is
     * number of hidden nodes = number of particles. */
    void setNumberOfHiddenUnits(const int numberOfHiddenUnits);

    /* Call this function from main.cpp to specify the number
     * of Metropolis cycles used in sampling. The default
     * setting is M=2^19=524,288. */
    void setNumberOfMetropolisCycles(const int steps);

    /* Call this function from main.cpp to specify the maximum
     * number of iterations used in the simulation. Default
     * setting is 1000. */
    void setNumberOfIterations(const int numberOfIterations = 1000);

    /* The adaptiveSteps function controls the number of cycles
     * used in each iteration. If adaptive steps is applied,
     * this function increases the number of cycles at the
     * right time. */
    int adaptiveSteps();

    /* The initializeSystem function initializes the system in the
     * right way. Various objects need to be initialized in the
     * correct order, which is the task of this function. This lets
     * the calls being in arbitrary order in main/configuration file. */
    void initializeSystem();

    /* The runMetropolisCycles function contains the Monte Carlo
     * sampling loop. This is the function can be trivially
     * parallelized. */
    void runMetropolisCycles();

    /* The checkingConvergence function is responsible for
     * checking whether the simulation has converged or not. */
    void checkingConvergence();

    /* The setNumberOfFreeDimensions function specifies the
     * number of free dimensions (the number of particles times
     * the number of dimensions). This is called automatically
     * when both the number of particles and number of dimensions
     * are determined. */
    void setNumberOfFreeDimensions();

    /* The setMaxParameters function speficies the maximum number
     * of parameters found in a single wave function element. This
     * is called automatically when all the wave function elements
     * are specified. */
    void setMaxParameters();

    /* The setGradients function declares the gradient matrix based
     * on the maximum number of parameters in a single wave function
     * element and the number of wave function elements. This is
     * called automatically when all the wave function elements are
     * specified. */
    void initializeMPI();
    void printLogo();
    void printInitialInformation();
    void printSystemInformation();
    void printHeaderLine();
    void setGradients();
    void setGlobalArraysToCalculate();
    void resetAllArrays();
    void setAllConstants();
    void collectAllLabels();
    void setAtomicNumber(const int Z);
    void runSimulation();
    void printToTerminal();
    void setNumberOfElements(const unsigned long numberOfElements);
    void initializeFromConfig(int argc, char** argv);
    void parser(const std::string configFile);
    void searchShortning(const std::vector<std::string> labels,
                         const std::string newLabel,
                         std::string &allLabels);

    void setStepLength(const double stepLength);
    void setEquilibrationFraction(const double equilibrationFraction);
    void setFrequency(const double omega);
    void setWidth(const double sigma);
    void setTotalSpin(const double totalSpin);
    void setLearningRate(const double eta);
    void setPath(const std::string path);

    void computeRadialOneBodyDensity(int numberOfBins = 1000, double maxRadius = 50);
    void computeSpatialOneBodyDensity(int numberOfBins = 1000, double maxRadius = 50);
    void computeTwoBodyDensity(int numberOfBins = 1000, double maxRadius = 50);
    void dumpEnergyToFile(bool printEnergyFile = true);
    void dumpParametersToFile(const bool printParametersToFile = true);
    void setConvergenceTools(const int numberOfEnergies = 5,
                             const double tolerance = 1e-6);
    void setAdaptiveStepTools(const int rangeOfAdaptiveSteps = 10,
                              const int additionalSteps = 4,
                              const int additionalStepsLastIteration = 8);
    void setScreeningTools(const double screeningStrength, const double dsl);
    void doBlocking(bool printInstantEnergyFile = true);

    void updateAllParameters(const Eigen::MatrixXd parameters);
    void initializeAllArrays(const Eigen::VectorXd positions,
                             const Eigen::VectorXd radialVector,
                             const Eigen::MatrixXd distanceMatrix);
    void updateAllArrays(const Eigen::VectorXd positions,
                         const Eigen::VectorXd radialVector,
                         const Eigen::MatrixXd distanceMatrix,
                         const int changedCoord);

    double evaluateProbabilityRatio();
    double getKineticEnergy();
    Eigen::MatrixXd getAllParameterGradients();

    void setHamiltonian(class Hamiltonian *hamiltonian);
    void setBasis(class Basis *basis);
    void setInitialState(class InitialState *initialState);
    void setInitialWeights(class InitialWeights *initialWeights);
    void setInteractionStyle(class Interaction *interaction);
    void setMetropolis(class Metropolis *metropolis);
    void setOptimization(class Optimization *optimization);
    void setRandomNumberGenerator(class RandomNumberGenerator *randomNumberGenerator);
    void setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements);
    void setWaveFunctionElement(WaveFunction *waveFunction);
    void setInputLayer(int numberOfUnits);
    void addDenseLayer(int numberOfUnits, Activation *activation);
    void setOutputLayer(Activation *activation);

    class WaveFunction *getWaveFunction() { return m_waveFunction; }
    class Hamiltonian *getHamiltonian() { return m_hamiltonian; }
    class Basis *getBasis() { return m_basis; }
    class Sampler *getSampler() { return m_sampler; }
    class Optimization *getOptimization() { return m_optimization; }
    class InitialWeights *getInitialWeights() { return m_initialWeights; }
    class InitialState *getInitialState() { return m_initialState; }
    class Interaction *getInteractionStyle() { return m_interactionStyle; }
    class RandomNumberGenerator *getRandomNumberGenerator() { return m_randomNumberGenerator; }

    int getNumberOfElements() { return m_numberOfElements; }
    int getNumberOfProcesses() { return m_numberOfProcesses; }
    int getNumberOfParticles() { return m_numberOfParticles; }
    int getNumberOfDimensions() { return m_numberOfDimensions; }
    int getNumberOfHiddenUnits() { return m_numberOfHiddenUnits; }
    int getNumberOfFreeDimensions() { return m_degreesOfFreedom; }
    int getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    int getMaxParameters() { return m_maxParameters; }
    int getAtomicNumber() { return m_Z; }
    int getRank() { return m_rank; }
    int getNumberOfBins() { return m_numberOfBins; }

    int getTotalStepsWOEqui() { return m_totalStepsWOEqui; }
    int getTotalStepsWEqui() { return m_totalStepsWEqui; }
    int getTotalEquilibriationSteps() { return m_totalEquilibriationSteps; }
    int getStepsWOEqui() { return m_stepsWOEqui; }
    int getStepsWEqui() { return m_stepsWEqui; }
    int getInitialTotalStepsWOEqui() { return m_initialTotalStepsWOEqui; }
    int getEquilibriationSteps() { return m_equilibriationSteps; }

    double getMaxRadius() { return m_maxRadius; }
    double getEquilibrationFraction() { return m_equilibrationFraction; }
    double getFrequency() { return m_omega; }
    double getWidth() { return m_sigma; }
    double getSorting() { return m_sorting; }
    double getLearningRate() { return m_eta; }
    double getStepLength() { return m_stepLength; }
    double getTotalSpin() { return m_totalSpin; }
    double getScreeningStrength() { return m_screeningStrength; }
    double getDSL() { return m_dsl; }

    bool getScreening() { return m_screening; }
    bool radialOneBodyDensity() { return m_computeOneBodyDensity; }
    bool spatialOneBodyDensity() { return m_computeOneBodyDensity2; }
    bool radialTwoBodyDensity() { return m_computeTwoBodyDensity; }
    bool printEnergyToFile() { return m_printEnergyToFile; }
    bool printParametersToFile() { return m_printParametersToFile; }
    bool doResampling() { return m_doResampling; }
    bool getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool getCalculateRadialVector() { return m_calculateRadialVector; }

    Eigen::VectorXd getPositions() { return m_positions; }
    Eigen::VectorXd getRadialVector() { return m_radialVector; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::MatrixXd getWeights() { return m_parameters; }

    std::string getPath() { return m_path; }
    std::string getTrialWaveFunction() { return m_trialWaveFunction; }

    std::vector<class WaveFunction *> getWaveFunctionElements() { return m_waveFunctionElements; }
    //std::vector<class Activation *> getActivationFunctions() { return m_activationFunctions; }
    std::vector<class Layer *> getLayers() { return m_layers; }
    //std::vector<int> getHiddenUnits() { return m_hiddenUnits; }

    // trim from left
    inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
    {
        s.erase(0, s.find_first_not_of(t));
        return s;
    }

    // trim from right
    inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
    {
        s.erase(s.find_last_not_of(t) + 1);
        return s;
    }

    // trim from left & right
    inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
    {
        return ltrim(rtrim(s, t), t);
    }

    // split string by whitespace
    std::vector<std::string> split(const std::string s) {
        std::istringstream iss(s);
        std::vector<std::string> splitted((std::istream_iterator<std::string>(iss)),
                                           std::istream_iterator<std::string>());
        return splitted;
    }

private:
    int m_numberOfParticles;
    int m_numberOfDimensions;
    int m_degreesOfFreedom;
    int m_numberOfIterations;
    int m_numberOfElements;
    int m_numberOfHiddenUnits;
    int m_maxParameters;
    int m_totalNumberOfParameters;
    int m_Z;
    int m_numberOfProcesses;
    int m_rank;

    int m_rangeOfAdaptiveSteps = 10;
    int m_additionalSteps = 4;
    int m_additionalStepsLastIter = 8;
    int m_numberOfNormalIterations = 1;
    int m_numberOfEnergies = 5;
    int m_numberOfBins = 1000;
    int m_iter = 0;
    int m_args = 1;
    int m_checkpointFreq = 100;

    int m_totalStepsWOEqui = int(pow(2,19));
    int m_totalStepsWEqui = int(pow(2,19));
    int m_totalEquilibriationSteps = 0;
    int m_stepsWOEqui = int(pow(2, 17));
    int m_stepsWEqui = int(pow(2, 17));
    int m_equilibriationSteps = 0;
    int m_initialStepsWOEqui = int(pow(2, 17));
    int m_initialTotalStepsWOEqui = int(pow(2, 19));

    double m_equilibrationFraction = 0;
    double m_tolerance = 1e-7;
    double m_maxRadius = 5;
    double m_totalTime = 0;
    double m_totalSpin = 0;
    double m_screeningStrength = 100;
    double m_dsl = 100;
    double m_globalTime = 0;

    double m_omega = 1.0;
    double m_sigma = 1.0;
    double m_stepLength = 0.05;
    double m_eta = 0.1;

    bool m_checkConvergence = false;
    bool m_applyAdaptiveSteps = false;
    bool m_computeOneBodyDensity = false;
    bool m_computeOneBodyDensity2 = false;
    bool m_computeTwoBodyDensity = false;
    bool m_printEnergyToFile = false;
    bool m_doResampling = false;
    bool m_printParametersToFile = false;
    bool m_calculateDistanceMatrix = false;
    bool m_calculateRadialVector = false;
    bool m_screening = false;
    bool m_sorting = true;

    std::chrono::system_clock::time_point m_start;

    class WaveFunction *m_waveFunction = nullptr;
    class Hamiltonian *m_hamiltonian = nullptr;
    class Basis *m_basis = new None(this);
    class InitialState *m_initialState = new RandomNormal(this);
    class InitialWeights *m_initialWeights = new Automatize(this);
    class Sampler *m_sampler = nullptr;
    class Metropolis *m_metropolis = new ImportanceSampling(this);
    class Optimization *m_optimization = new ADAM(this);
    class RandomNumberGenerator *m_randomNumberGenerator = new MersenneTwister();
    class Interaction *m_interactionStyle = new NoInteraction(this);

    std::vector<class WaveFunction *> m_waveFunctionElements;
    std::vector<class Layer *> m_layers;
    //std::vector<class Activation *> m_activationFunctions;
    //std::vector<int> m_hiddenUnits;

    std::string m_path = "";
    std::string m_trialWaveFunction;
    std::string m_configFile = "none";

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_gradients;
    Eigen::VectorXd m_energies;
};
