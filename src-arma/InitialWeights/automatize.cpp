#include "automatize.h"
#include "../Optimization/optimization.h"
#include "../system.h"
#include "constant.h"
#include "randomize.h"
#include <fstream>
#include <iostream>

Automatize::Automatize(System *system)
    : InitialWeights(system)
{
    m_system = system;
}

std::string Automatize::generateFileName(std::string name, std::string extension)
{
    std::string fileName = m_path;
    /*
    fileName += "arma::uword" + std::to_string(m_interaction) + "/";
    fileName += m_hamiltonian + "/";
    fileName += name + "/";
    fileName += m_trialWaveFunction + "/";
    fileName += std::to_string(m_numberOfDimensions) + "D/";
    fileName += std::to_string(m_numberOfParticles) + "P/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += m_system->getOptimization()->getLabel();
    fileName += "_MC" + std::to_string(m_initialTotalStepsWOEqui);
    fileName += extension;
    */
    return fileName;
}

void writeFileContent(std::ifstream infile, arma::mat &matrix)
{
    std::string line;
    double value;
    char delimiter;
    arma::uword i = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        arma::uword j = 0;
        while (iss >> value) {
            matrix(i, j) = value;
            iss >> delimiter;
            j++;
        }
        i++;
    }
}

void Automatize::setupInitialWeights()
{
    /*
    bool searchForWeights = true;
    std::ifstream infile(generateFileName("weights", ".dat"));
    if (infile.is_open() && searchForWeights) {
        m_parameters = arma::mat::Zero(m_numberOfElements, m_maxParameters);
        std::string line;
        double value;
        char delimiter;
        arma::uword i = 0;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            arma::uword j = 0;
            while (iss >> value) {
                m_parameters(i, j) = value;
                iss >> delimiter;
                j++;
            }
            i++;
        }
    } else {
    */

    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_trialWaveFunction = m_system->getTrialWaveFunction();

    if (m_trialWaveFunction == "VMC") {
        m_method = new Constant(m_system, 1.0);
    } else if (m_trialWaveFunction == "RBM") {
        m_method = new Randomize(m_system, 0.2);
    } else if (m_trialWaveFunction == "RBMSJ") {
        m_method = new Randomize(m_system, 0.2);
    } else if (m_trialWaveFunction == "RBMPJ") {
        m_method = new Constant(m_system, 0.0);
    } else if (m_trialWaveFunction == "PRBM") {
        m_method = new Constant(m_system, 0.0);
    } else {
        m_method = new Randomize(m_system, 0.01);
    }
    m_method->setupInitialWeights();
    m_parameters = m_method->getParameters();
    //}
    m_system->updateAllParameters(m_parameters);
}

arma::mat Automatize::getParameters()
{
    return m_parameters;
}