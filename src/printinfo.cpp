#include <iostream>
#include <string>
#include <ctime>
#include <chrono>


#include "system.h"


/* ----------------------------------------------------------------------------
  Print logo to terminal
---------------------------------------------------------------------------- */

void System::printLogo()
{
    std::cout << "__      ____  __          _____ _    _ _____ _   _ ______ " << std::endl;
    std::cout << "\\ \\    / /  \\/  |   /\\   / ____| |  | |_   _| \\ | |  ____|" << std::endl;
    std::cout << " \\ \\  / /| \\  / |  /  \\ | |    | |__| | | | |  \\| | |__   " << std::endl;
    std::cout << "  \\ \\/ / | |\\/| | / /\\ \\| |    |  __  | | | | . ` |  __|  " << std::endl;
    std::cout << "   \\  /  | |  | |/ ____ \\ |____| |  | |_| |_| |\\  | |____ " << std::endl;
    std::cout << "    \\/   |_|  |_/_/    \\_\\_____|_|  |_|_____|_| \\_|______|" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
  Print initial information to terminal
---------------------------------------------------------------------------- */

void System::printInitialInformation()
{
    m_start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(m_start);
    std::cout << "Started computation at " << std::ctime(&start_time)
              << "Running on " << m_numberOfProcesses << " CPU threads using Open MPI" << std::endl;
    std::cout << "Simulation is run from the directory: " << m_path << std::endl;
}


/* ----------------------------------------------------------------------------
  Print system information to terminal
---------------------------------------------------------------------------- */

void System::printSystemInformation()
{
    std::cout << std::fixed;
    std::cout << std::boolalpha;
    std::cout << std::setprecision(6);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "             SYSTEM INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of particles:      " << m_numberOfParticles << std::endl;
    std::cout << "Number of dimensions:     " << m_numberOfDimensions << std::endl;
    std::cout << "Interaction style:        " << m_interactionStyle->getLabel() << std::endl;
    std::cout << "Hamiltonian:              " << m_hamiltonian->getLabel() << std::endl;
    std::cout << "Initial state:            " << m_initialState->getLabel() << std::endl;
    if (m_hamiltonian->getLabel() == "harmonic oscillator") {
        std::cout << "Oscillator frequency:     " << m_omega << std::endl;
    } else if (m_hamiltonian->getLabel() == "double well"){
        std::cout << "Oscillator frequency:     " << m_omega << std::endl;
    } else if (m_hamiltonian->getLabel() == "atom") {
        std::cout << "Atomic number:            " << m_Z << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "           WAVE FUNCTION INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    for (int i = 0; i < m_numberOfElements; i++) {
        std::cout << "Element " << i << ":                " << m_waveFunctionElements[unsigned(i)]->getLabel() << std::endl;
    }
    std::cout << "Basis:                    " << m_basis->getLabel() << std::endl;
    std::cout << "Initial parameters:       " << m_initialWeights->getLabel() << std::endl;
    std::cout << "Number of parameters:     " << this->getTotalNumberOfParameters() << std::endl;
    std::cout << "Number of hidden nodes:   " << m_numberOfHiddenUnits << std::endl;
    std::cout << "Print parameters to file: " << m_printParametersToFile << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "            SIMULATION INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Max number of iterations: " << m_numberOfIterations << std::endl;
    std::cout << "Number of MC cycles:      " << m_totalStepsWOEqui << std::endl;
    std::cout << "Number of burn-in cycles: " << m_burnInSteps << std::endl;
    std::cout << "Learning rate:            " << m_eta << std::endl;
    std::cout << "Step length:              " << m_stepLength << std::endl;
    //std::cout << "Equilibration fraction:   " << m_equilibrationFraction << std::endl;
    std::cout << "Sampling:                 " << m_metropolis->getLabel() << std::endl;
    std::cout << "Optimization:             " << m_optimization->getLabel() << std::endl;
    std::cout << "Random number generator:  " << m_randomNumberGenerator->getLabel() << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "               PARTICLE DENSITY" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Compute radial one-body density:  " << m_computeOneBodyDensity << std::endl;
    std::cout << "Compute spatial one-body density: " << m_computeOneBodyDensity2 << std::endl;
    std::cout << "Compute radial two-body density:  " << m_computeTwoBodyDensity << std::endl;
    std::cout << "Max radius of electron density:   " << m_maxRadius << std::endl;
    std::cout << "Number of bins in each direction: " << m_numberOfBins << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                  RESAMPLING" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Do resampling:            " << m_doResampling << std::endl;
    std::cout << "Print energies to file:   " << m_printEnergyToFile << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                  CONVERGENCE" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Check convergence:        " << m_checkConvergence << std::endl;
    std::cout << "Tolerance:                " << m_tolerance << std::endl;
    std::cout << "Number of energies:       " << m_numberOfEnergies << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                 ADAPTIVE STEPS" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Apply adaptive steps:       " << m_applyAdaptiveSteps << std::endl;
    std::cout << "Range of adaptive steps:    " << m_rangeOfAdaptiveSteps << std::endl;
    std::cout << "Additional steps:           " << m_additionalSteps << std::endl;
    std::cout << "Additional steps last iter: " << m_additionalStepsLastIter << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
  Print header line to terminal, containing information about the various
  columns
---------------------------------------------------------------------------- */

void System::printHeaderLine()
{
    std::cout << "Step  " << "Energy  " << "Energy_STD  " << "Kinetic  "
              << "Kinetic_STD  " << "External  " << "External_STD  "
              << "Interaction  " << "Interaction_STD  " << "Acceptence  "
              << "CPU_time" << std::endl;
}

