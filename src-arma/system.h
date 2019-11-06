#pragma once
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <armadillo>

#include "main.h"

class System
{
public:
    arma::uword adaptiveSteps();
    void initializeSystem();
    void runMetropolisCycles();
    void checkingConvergence();
    void setNumberOfFreeDimensions();
    void setMaxParameters();
    void setGradients();
    void setGlobalArraysToCalculate();
    void resetAllArrays();
    void setAllConstants();
    void collectAllLabels();
    void setAtomicNumber(const arma::uword Z);
    void runSimulation(const arma::uword numberOfIterations = 1000);
    void printToTerminal(const arma::uword numberOfIterations);
    void setNumberOfParticles(const arma::uword numberOfParticles);
    void setNumberOfDimensions(const arma::uword numberOfDimensions);
    void setNumberOfHiddenNodes(const arma::uword numberOfHiddenNodes);
    void setNumberOfMetropolisSteps(const arma::uword steps);
    void setNumberOfElements(const unsigned long numberOfElements);
    void setStepLength(const double stepLength);
    void setEquilibrationFraction(const double equilibrationFraction);
    void setFrequency(const double omega);
    void setWidth(const double sigma);
    void setTotalSpin(const double totalSpin);
    void setLearningRate(const double eta);
    void setPath(const std::string path);
    void parserConstants(const std::string configFile, arma::uword &numberOfIterations);
    void parserObjects(const std::string configFile);
    void searchShortning(const std::vector<std::string> labels,
                         const std::string newLabel,
                         std::string &allLabels);

    void setInteraction(const bool interaction);
    void setParameterPrintingTools(const bool printParametersToFile);
    void setConvergenceTools(const bool checkConvergence,
                             const arma::uword numberOfEnergies,
                             const double tolerance);
    void setAdaptiveStepTools(const bool applyAdaptiveSteps,
                              const arma::uword rangeOfAdaptiveSteps,
                              const arma::uword additionalSteps,
                              const arma::uword additionalStepsLastIteration);
    void setDensityTools(const bool computeOneBodyDensity,
                         const bool computeOneBodyDensity2,
                         const bool computeTwoBodyDensity,
                         const arma::uword numberOfBins,
                         const double maxRadius);
    void setScreeningTools(const bool screening, const double screeningStrength, const double dsl);
    void setEnergyPrintingTools(const bool printEnergyFile, const bool printInstantEnergyFile);
    void initializeMPI();

    void updateAllParameters(const arma::mat parameters);
    void initializeAllArrays(const arma::vec positions,
                             const arma::vec radialVector,
                             const arma::mat distanceMatrix);
    void updateAllArrays(const arma::vec positions,
                         const arma::vec radialVector,
                         const arma::mat distanceMatrix,
                         const arma::uword changedCoord);

    double evaluateProbabilityRatio();
    double getKineticEnergy();
    arma::mat getAllParameterGradients();

    void setHamiltonian(class Hamiltonian *hamiltonian);
    void setBasis(class Basis *basis);
    void setInitialState(class InitialState *initialState);
    void setInitialWeights(class InitialWeights *initialWeights);
    void setMetropolis(class Metropolis *metropolis);
    void setOptimization(class Optimization *optimization);
    void setRandomNumberGenerator(class RandomNumberGenerator *randomNumberGenerator);
    void setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements);
    void setWaveFunctionElement(WaveFunction *waveFunction);

    class WaveFunction *getWaveFunction() { return m_waveFunction; }
    class Hamiltonian *getHamiltonian() { return m_hamiltonian; }
    class Basis *getBasis() { return m_basis; }
    class Sampler *getSampler() { return m_sampler; }
    class Optimization *getOptimization() { return m_optimization; }
    class InitialWeights *getInitialWeights() { return m_initialWeights; }
    class InitialState *getInitialState() { return m_initialState; }
    class RandomNumberGenerator *getRandomNumberGenerator() { return m_randomNumberGenerator; }

    arma::uword getNumberOfElements() { return m_numberOfElements; }
    arma::uword getNumberOfProcesses() { return m_numberOfProcesses; }
    arma::uword getNumberOfParticles() { return m_numberOfParticles; }
    arma::uword getNumberOfDimensions() { return m_numberOfDimensions; }
    arma::uword getNumberOfHiddenNodes() { return m_numberOfHiddenNodes; }
    arma::uword getNumberOfFreeDimensions() { return m_degreesOfFreedom; }
    arma::uword getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    arma::uword getMaxParameters() { return m_maxParameters; }
    arma::uword getAtomicNumber() { return m_Z; }
    arma::uword getRank() { return m_myRank; }
    arma::uword getNumberOfBins() { return m_numberOfBins; }

    arma::uword getTotalStepsWOEqui() { return m_totalStepsWOEqui; }
    arma::uword getTotalStepsWEqui() { return m_totalStepsWEqui; }
    arma::uword getTotalEquilibriationSteps() { return m_totalEquilibriationSteps; }
    arma::uword getStepsWOEqui() { return m_stepsWOEqui; }
    arma::uword getStepsWEqui() { return m_stepsWEqui; }
    arma::uword getInitialTotalStepsWOEqui() { return m_initialTotalStepsWOEqui; }
    arma::uword getEquilibriationSteps() { return m_equilibriationSteps; }

    double getMaxRadius() { return m_maxRadius; }
    double getEquilibrationFraction() { return m_equilibrationFraction; }
    double getFrequency() { return m_omega; }
    double getWidth() { return m_sigma; }
    double getLearningRate() { return m_eta; }
    double getStepLength() { return m_stepLength; }
    double getTotalSpin() { return m_totalSpin; }
    double getScreeningStrength() { return m_screeningStrength; }
    double getDSL() { return m_dsl; }
    bool getScreening() { return m_screening; }
    bool getInteraction() { return m_interaction; }
    bool computeOneBodyDensity() { return m_computeOneBodyDensity; }
    bool computeOneBodyDensity2() { return m_computeOneBodyDensity2; }
    bool computeTwoBodyDensity() { return m_computeTwoBodyDensity; }
    bool printEnergyToFile() { return m_printEnergyToFile; }
    bool printParametersToFile() { return m_printParametersToFile; }
    bool doResampling() { return m_doResampling; }
    bool getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool getCalculateRadialVector() { return m_calculateRadialVector; }
    arma::vec getPositions() { return m_positions; }
    arma::vec getRadialVector() { return m_radialVector; }
    arma::mat getDistanceMatrix() { return m_distanceMatrix; }
    arma::mat getWeights() { return m_parameters; }
    std::string getPath() { return m_path; }
    std::string getTrialWaveFunction() { return m_trialWaveFunction; }
    std::vector<class WaveFunction *> getWaveFunctionElements() { return m_waveFunctionElements; }

private:
    arma::uword m_numberOfParticles;
    arma::uword m_numberOfDimensions;
    arma::uword m_degreesOfFreedom;
    arma::uword m_numberOfElements;
    arma::uword m_numberOfHiddenNodes;
    arma::uword m_maxParameters;
    arma::uword m_totalNumberOfParameters;
    arma::uword m_Z;
    int m_numberOfProcesses;
    int m_myRank;

    arma::uword m_rangeOfAdaptiveSteps = 10;
    arma::uword m_additionalSteps = 4;
    arma::uword m_additionalStepsLastIter = 8;
    arma::uword m_lastIteration = 1;
    arma::uword m_numberOfEnergies = 5;
    arma::uword m_numberOfBins = 1000;
    arma::uword m_iter = 0;

    arma::uword m_totalStepsWOEqui = arma::uword(pow(2,19));
    arma::uword m_totalStepsWEqui = arma::uword(pow(2,19));
    arma::uword m_totalEquilibriationSteps = 0;
    arma::uword m_stepsWOEqui = arma::uword(pow(2, 17));
    arma::uword m_stepsWEqui = arma::uword(pow(2, 17));
    arma::uword m_equilibriationSteps = 0;
    arma::uword m_initialStepsWOEqui = arma::uword(pow(2, 17));
    arma::uword m_initialTotalStepsWOEqui = arma::uword(pow(2, 19));

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

    bool m_interaction = true;
    bool m_checkConvergence = false;
    bool m_applyAdaptiveSteps = false;
    bool m_computeOneBodyDensity = false;
    bool m_computeOneBodyDensity2 = false;
    bool m_computeTwoBodyDensity = false;
    bool m_printEnergyToFile = false;
    bool m_doResampling = true;
    bool m_printParametersToFile = false;
    bool m_calculateDistanceMatrix = false;
    bool m_calculateRadialVector = false;
    bool m_screening = false;

    class WaveFunction *m_waveFunction = nullptr;
    class Hamiltonian *m_hamiltonian = new HarmonicOscillator(this);
    class Basis *m_basis = new Hermite(this);
    class InitialState *m_initialState = new RandomNormal(this);
    class InitialWeights *m_initialWeights = new Automatize(this);
    class Sampler *m_sampler = nullptr;
    class Metropolis *m_metropolis = new ImportanceSampling(this);
    class Optimization *m_optimization = new ADAM(this);
    class RandomNumberGenerator *m_randomNumberGenerator = new MersenneTwister();
    std::vector<class WaveFunction *> m_waveFunctionElements;

    std::string m_path = "../data/";
    std::string m_trialWaveFunction;

    arma::vec m_positions;
    arma::mat m_parameters;
    arma::mat m_distanceMatrix;
    arma::vec m_radialVector;
    arma::mat m_gradients;
    arma::vec m_energies;
};
