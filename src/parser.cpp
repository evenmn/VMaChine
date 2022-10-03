#include <iostream>
#include <sstream>
#include <fstream>
#include <Eigen/Dense>

#include "system.h"


/* ----------------------------------------------------------------------------
  Initialize variables from configuration file. Will overwrite main and
  default settings
---------------------------------------------------------------------------- */

void System::initializeFromConfig(int argc, char** argv) {
    if (argc >= 2) {
        m_configFile = argv[1];
        m_args = argc;
    }
}

/* ----------------------------------------------------------------------------
  Parse input script
---------------------------------------------------------------------------- */

void System::parser(const std::string configFile)
{
    std::ifstream infile;
    infile.open(configFile.c_str());
    if (!infile.is_open() && m_args >= 2) {
        std::cout << std::endl;
        std::cerr << "File: '" << configFile << "'" << std::endl;
        perror("File not found");
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
    else {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream is_line(line);
            std::string key;
            if (line.empty()) {
                /* Continue if line is blank */
                continue;
            } else if (line.rfind("#", 0) == 0) {
                /* Continue if line starts with '#' */
                continue;
            } else if (std::getline(is_line, key, ':')) {
                key = trim(key);
                std::string value;
                if (std::getline(is_line, value)) {
                    std::vector<std::string> splitted = split(value);
                    if (key == "numParticles") {
                        m_numberOfParticles = std::stoi(splitted.at(0));
                        m_numberOfHiddenUnits = m_numberOfParticles;
                        m_Z = m_numberOfParticles;
                    } else if (key == "numDimensions") {
                        m_numberOfDimensions = std::stoi(splitted.at(0));
                    } else if (key == "omega") {
                        m_omega = std::stod(splitted.at(0));
                        m_stepLength = 0.1 / sqrt(m_omega);
                        m_sigma = 1.0 / sqrt(m_omega);
                    } else if (key == "atomicNumber") {
                        m_Z = std::stoi(splitted.at(0));
                    } else if (key == "learningRate") {
                        m_eta = std::stod(splitted.at(0));
                    } else if (key == "maxRadius") {
                        m_maxRadius = std::stod(splitted.at(0));
                    } else if (key == "numIterations") {
                        m_numberOfIterations = std::stoi(splitted.at(0));
                    //} else if (key == "equilibration") {
                    //    setEquilibrationFraction(std::stod(splitted.at(0)));
                    } else if (key == "burnin") {
                        setBurnInSteps(std::stoi(splitted.at(0)));
                    } else if (key == "numSteps") {
                        setNumberOfMetropolisCycles(std::stoi(splitted.at(0)));
                    } else if (key == "numHiddenNodes") {
                        m_numberOfHiddenUnits = std::stoi(splitted.at(0));
                    } else if (key == "sorting") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_sorting;
                    } else if (key == "totalSpin") {
                        m_totalSpin = std::stod(splitted.at(0));
                    } else if (key == "stepLength") {
                        m_stepLength = std::stod(splitted.at(0));
                    } else if (key == "checkpointFreq") {
                        m_checkpointFreq = std::stod(splitted.at(0));
                    } else if (key == "checkConvergence") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_checkConvergence;
                    } else if (key == "applyAdaptiveSteps") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_applyAdaptiveSteps;
                    } else if (key == "computeOneBodyDensity") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeOneBodyDensity;
                    } else if (key == "computeOneBodyDensity2") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeOneBodyDensity2;
                    } else if (key == "computeTwoBodyDensity") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeTwoBodyDensity;
                    } else if (key == "printEnergyToFile") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_printEnergyToFile;
                    } else if (key == "printParametersToFile") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_printParametersToFile;
                    } else if (key == "doResampling") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_doResampling;
                    } else if (key == "path") {
                        m_path = splitted.at(0);
                    } else if (key == "numberOfEnergies") {
                        m_numberOfEnergies = std::stoi(splitted.at(0));
                    } else if (key == "tolerance") {
                        m_tolerance = std::stod(splitted.at(0));
                    } else if (key == "rangeOfAdaptiveSteps") {
                        m_rangeOfAdaptiveSteps = std::stoi(splitted.at(0));
                    } else if (key == "additionalSteps") {
                        m_additionalSteps = std::stoi(splitted.at(0));
                    } else if (key == "additionalStepsLastIter") {
                        m_additionalStepsLastIter = std::stoi(splitted.at(0));
                    } else if (key == "numberOfBins") {
                        m_numberOfBins = std::stoi(splitted.at(0));
                    } else if (key == "basis") {
                        delete m_basis;
                        if (splitted.at(0) == "hermite") {
                            setBasis(new Hermite(this));
                        //} else if (splitted.at(0) == "hermiteExpansion") {
                        //    setBasis(new HermiteExpansion(this));
                        } else if (splitted.at(0) == "hydrogenOrbital") {
                            setBasis(new HydrogenOrbital(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Basis '"
                                      << splitted.at(0)
                                      << "' is not implemented!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "hamiltonian") {
                        delete m_hamiltonian;
                        if (splitted.at(0) == "harmonicOscillator") {
                            setHamiltonian(new HarmonicOscillator(this));
                        } else if (splitted.at(0) == "doubleWell") {
                            setHamiltonian(new DoubleWell(this, std::stod(splitted.at(1))));
                        } else if (splitted.at(0) == "atomicNucleus") {
                            setHamiltonian(new AtomicNucleus(this));
                        } else if (splitted.at(0) == "ellipticalHarmonicOscillator") {
                            setHamiltonian(new EllipticalHarmonicOscillator(this, std::stod(splitted.at(1))));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "The Hamiltonian '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "optimization") {
                        delete m_optimization;
                        if (splitted.at(0) == "adam") {
                            setOptimization(new ADAM(this));
                        } else if (splitted.at(0) == "gd") {
                            if (splitted.size() >= 3) {
                                setOptimization(new GradientDescent(this, std::stod(splitted.at(1)), std::stod(splitted.at(2))));
                            } else {
                                setOptimization(new GradientDescent(this, 0.0, 0.0));
                            }
                        } else if (splitted.at(0) == "sgd") {
                          if (splitted.size() >= 3) {
                              setOptimization(new SGD(this, std::stod(splitted.at(1)), std::stod(splitted.at(2))));
                          } else {
                              setOptimization(new SGD(this, 0.0, 0.0));
                          }
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Optimization method '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "initialWeights") {
                        delete m_initialWeights;
                        if (splitted.at(0) == "automatize") {
                            setInitialWeights(new Automatize(this));
                        } else if (splitted.at(0) == "randomuniform") {
                            if (splitted.size() >= 2) {
                                setInitialWeights(new RandomUniformWeights(this, std::stod(splitted.at(1))));
                            } else {
                                setInitialWeights(new RandomUniformWeights(this));
                            }
                        } else if (splitted.at(0) == "constant") {
                            if (splitted.size() >= 2) {
                                setInitialWeights(new Constant(this, std::stod(splitted.at(1))));
                            } else {
                                setInitialWeights(new Constant(this, 1.0));
                            }
                        } else if (splitted.at(0) == "fromfile") {
                            if (splitted.size() >= 2) {
                                setInitialWeights(new FromFile(this, splitted.at(1)));
                            } else {
                                setInitialWeights(new FromFile(this));
                            }
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Initial parameter configuration '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "initialState") {
                        delete m_initialState;
                        if (splitted.at(0) == "randomNormal") {
                            setInitialState(new RandomNormal(this));
                        } else if (splitted.at(0) == "randomUniform") {
                            setInitialState(new RandomUniform(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Initial state configuration '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "sampling") {
                        delete m_metropolis;
                        if (splitted.at(0) == "importanceSampling") {
                            setMetropolis(new ImportanceSampling(this));
                        } else if (splitted.at(0) == "bruteForce") {
                            setMetropolis(new BruteForce(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Sampling method '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "waveFunction") {
                        for (int i=m_numberOfElements; i--;)
                        {
                            delete m_waveFunctionElements[i];
                        }
                        m_waveFunctionElements.clear();

                        std::vector<WaveFunction *> waveFunctionElements;
                        if (splitted.at(0) == "VMC") {
                            m_waveFunctionElements.push_back(new class Gaussian(this));
                            //waveFunctionElements.push_back(new class Gaussian(this));
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "RBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                        } else if (splitted.at(0) == "RBMPJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "RBMSJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class SimpleJastrow(this));
                        } else if (splitted.at(0) == "PRBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PartlyRestricted(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Wave function '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                        setWaveFunctionElements(waveFunctionElements);
                    } else if (key == "waveFunctionElement") {
                        if (splitted.at(0) == "gaussian") {
                            setWaveFunctionElement(new class Gaussian(this));
                        } else if (splitted.at(0) == "slaterDeterminant") {
                            setWaveFunctionElement(new class SlaterDeterminant(this));
                        } else if (splitted.at(0) == "padeJastrow") {
                            setWaveFunctionElement(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "simpleJastrow") {
                            setWaveFunctionElement(new class SimpleJastrow(this));
                        } else if (splitted.at(0) == "RBMGaussian") {
                            setWaveFunctionElement(new class RBMGaussian(this));
                        } else if (splitted.at(0) == "RBMProduct") {
                            setWaveFunctionElement(new class RBMProduct(this));
                        } else if (splitted.at(0) == "hydrogenLike") {
                            setWaveFunctionElement(new class HydrogenLike(this));
                        } else if (splitted.at(0) == "hardCoreJastrow") {
                            setWaveFunctionElement(new class HardCoreJastrow(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Wave function element '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "interactionStyle") {
                        delete m_interactionStyle;
                        if (splitted.at(0) == "noInteraction") {
                            setInteractionStyle(new class NoInteraction(this));
                        } else if (splitted.at(0) == "coulomb") {
                            setInteractionStyle(new class Coulomb(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Interaction style '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "rng") {
                        delete m_randomNumberGenerator;
                        if (splitted.at(0) == "MersenneTwister") {
                            setRandomNumberGenerator(new class MersenneTwister());
                            if (splitted.size() > 1) {
                                m_randomNumberGenerator->setSeed(std::stoi(splitted.at(1)));
                            }
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Random number generator '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else {
                        std::cout << std::endl;
                        std::cerr << "Invalid key '"
                                  << key
                                  << "' is passed to configuration file!" << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 143);
                    }
                }
            } else {
                std::cout << std::endl;
                std::cerr << "Invalid object detected in configuration file!" << std::endl;
                std::cerr << "Error raised when tried to read: "
                          << line
                          << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 143);
            }
        }
    }
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}


/* ----------------------------------------------------------------------------
  Map a wave function abbreviation to the actual wave function elements
---------------------------------------------------------------------------- */

void System::searchShortning(const std::vector<std::string> labels,
                             const std::string newLabel,
                             std::string &allLabels)
{
    if (labels.size() == 1) {
        if (allLabels == "_" + labels.at(0)) {
            allLabels = newLabel;
        }
    }
    if (labels.size() == 2) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                if (i != j) {
                    if (allLabels == "_" + i + "_" + j) {
                        allLabels = newLabel;
                        break;
                    }
                }
            }
        }
    } else if (labels.size() == 3) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                for (auto &k : labels) {
                    if (i != j || i != k || j != k) {
                        if (allLabels == "_" + i + "_" + j + "_" + k) {
                            allLabels = newLabel;
                            break;
                        }
                    }
                }
            }
        }
    } else if (labels.size() == 4) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                for (auto &k : labels) {
                    for (auto &l : labels) {
                        if (i != j || i != k || i != l || j != k || j != l || k != l) {
                            if (allLabels == "_" + i + "_" + j + "_" + k + "_" + l) {
                                allLabels = newLabel;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}


/* ----------------------------------------------------------------------------
  An overview of the possible trial wave function abbreviations
---------------------------------------------------------------------------- */

void System::collectAllLabels()
{
    m_trialWaveFunction = "";
    for (auto &i : m_waveFunctionElements) {
        m_trialWaveFunction += "_";
        m_trialWaveFunction += i->getLabel();
    }

    std::vector<std::string> testVMC1;
    testVMC1.push_back("Gaussian");
    searchShortning(testVMC1, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC2;
    testVMC2.push_back("Gaussian");
    testVMC2.push_back("Padé-Jastrow");
    searchShortning(testVMC2, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC3;
    testVMC3.push_back("Gaussian");
    testVMC3.push_back("Padé-Jastrow");
    testVMC3.push_back("Slater determinant");
    searchShortning(testVMC3, "VMC", m_trialWaveFunction);

    std::vector<std::string> testRBM1;
    testRBM1.push_back("RBM-Gaussian");
    testRBM1.push_back("RBM-product");
    searchShortning(testRBM1, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBM2;
    testRBM2.push_back("RBM-Gaussian");
    testRBM2.push_back("RBM-product");
    testRBM2.push_back("Slater determinant");
    searchShortning(testRBM2, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ1;
    testRBMPJ1.push_back("RBM-Gaussian");
    testRBMPJ1.push_back("RBM-product");
    testRBMPJ1.push_back("Padé-Jastrow");
    searchShortning(testRBMPJ1, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ2;
    testRBMPJ2.push_back("RBM-Gaussian");
    testRBMPJ2.push_back("RBM-product");
    testRBMPJ2.push_back("Padé-Jastrow");
    testRBMPJ2.push_back("Slater determinant");
    searchShortning(testRBMPJ2, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testPRBM1;
    testPRBM1.push_back("RBM-Gaussian");
    testPRBM1.push_back("RBM-product");
    testPRBM1.push_back("partlyrestricted");
    searchShortning(testPRBM1, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testPRBM2;
    testPRBM2.push_back("RBM-Gaussian");
    testPRBM2.push_back("RBM-product");
    testPRBM2.push_back("partlyrestricted");
    testPRBM2.push_back("Slater determinant");
    searchShortning(testPRBM2, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM1;
    testDRBM1.push_back("RBM-Gaussian");
    testDRBM1.push_back("drbmproduct");
    searchShortning(testDRBM1, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM2;
    testDRBM2.push_back("RBM-Gaussian");
    testDRBM2.push_back("drbmproduct");
    testDRBM2.push_back("Slater determinant");
    searchShortning(testDRBM2, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ1;
    testRBMSJ1.push_back("RBM-Gaussian");
    testRBMSJ1.push_back("RBM-product");
    testRBMSJ1.push_back("simple Jastrow");
    searchShortning(testRBMSJ1, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ2;
    testRBMSJ2.push_back("RBM-Gaussian");
    testRBMSJ2.push_back("RBM-product");
    testRBMSJ2.push_back("simple Jastrow");
    testRBMSJ2.push_back("Slater determinant");
    searchShortning(testRBMSJ2, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testVMC4;
    testVMC4.push_back("hydrogenlike");
    searchShortning(testVMC4, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC5;
    testVMC5.push_back("Slater determinant");
    searchShortning(testVMC5, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC6;
    testVMC6.push_back("hydrogenlike");
    testVMC6.push_back("Padé-Jastrow");
    searchShortning(testVMC6, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC7;
    testVMC7.push_back("Slater determinant");
    testVMC7.push_back("Padé-Jastrow");
    searchShortning(testVMC7, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC8;
    testVMC8.push_back("Slater determinant");
    testVMC8.push_back("Gaussian");
    searchShortning(testVMC8, "VMC", m_trialWaveFunction);

    std::vector<std::string> testSSJ1;
    testSSJ1.push_back("Slater determinant");
    testSSJ1.push_back("Gaussian");
    testSSJ1.push_back("simple Jastrow");
    searchShortning(testSSJ1, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testSSJ2;
    testSSJ2.push_back("Gaussian");
    testSSJ2.push_back("simple Jastrow");
    searchShortning(testSSJ2, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testFNN;
    testFNN.push_back("fnn");
    searchShortning(testFNN, "FNN", m_trialWaveFunction);

    std::vector<std::string> testBVMC;
    testBVMC.push_back("Gaussian");
    testBVMC.push_back("hard-core Jastrow");
    searchShortning(testBVMC, "BVMC", m_trialWaveFunction);
}
