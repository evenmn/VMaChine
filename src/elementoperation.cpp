#include <iostream>
#include <Eigen/Dense>

#include "system.h"

/* ----------------------------------------------------------------------------
  Initialize the wave function elements with essential variables specified by
  user
---------------------------------------------------------------------------- */

void System::setAllConstants()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_waveFunctionElements[unsigned(i)]->setConstants(i);
    }
}


/* ----------------------------------------------------------------------------
  Initialize the wave function elements with position
---------------------------------------------------------------------------- */

void System::initializeAllArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd radialVector,
                                 const Eigen::MatrixXd distanceMatrix)
{
    for (auto &i : m_waveFunctionElements) {
        i->initializeArrays(positions, radialVector, distanceMatrix);
        i->setArrays();
    }
}


/* ----------------------------------------------------------------------------
  Update positions in all wave function elements when a particle is moved
---------------------------------------------------------------------------- */

void System::updateAllArrays(const Eigen::VectorXd positions,
                             const Eigen::VectorXd radialVector,
                             const Eigen::MatrixXd distanceMatrix,
                             const int changedCoord)
{
    for (auto &i : m_waveFunctionElements) {
        i->setArrays();
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}


/* ----------------------------------------------------------------------------
  Reset positions in all wave function elements when a particle move is
  rejected
---------------------------------------------------------------------------- */

void System::resetAllArrays()
{
    for (auto &i : m_waveFunctionElements) {
        i->resetArrays();
    }
}


/* ----------------------------------------------------------------------------
  Update the parameters/weights in all the wave function element when the
  parameters are updated
---------------------------------------------------------------------------- */

void System::updateAllParameters(const Eigen::MatrixXd parameters)
{
    for (auto &i : m_waveFunctionElements) {
        i->updateParameters(parameters);
    }
}


/* ----------------------------------------------------------------------------
  Evaluate the collective probability ratio
---------------------------------------------------------------------------- */

double System::evaluateProbabilityRatio()
{
    double ratio = 1;
    for (auto &i : m_waveFunctionElements) {
        ratio *= i->evaluateRatio();
    }
    return ratio;
}


/* ----------------------------------------------------------------------------
  Obtain the total kinetic energy of the system based on the gradients and 
  Laplacians of all elements
---------------------------------------------------------------------------- */

double System::getKineticEnergy()
{
    double kineticEnergy = 0;
    for (auto &i : m_waveFunctionElements) {
        kineticEnergy += i->computeLaplacian();
    }
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        double nablaLnPsi = 0;
        for (auto &i : m_waveFunctionElements) {
            nablaLnPsi += i->computeGradient(k);
        }
        kineticEnergy += nablaLnPsi * nablaLnPsi;
    }
    return -0.5 * kineticEnergy;
}


/* ----------------------------------------------------------------------------
  Get the gradient of all the wave function elements with respect to all the 
  parameters. To be used in the parameter update
---------------------------------------------------------------------------- */

Eigen::MatrixXd System::getAllParameterGradients()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_gradients.row(i) = m_waveFunctionElements[unsigned(i)]->computeParameterGradient();
    }
    return m_gradients;
}


/* ----------------------------------------------------------------------------
  Check if the elements need distance matrix and/or radial distance vector
---------------------------------------------------------------------------- */

void System::setGlobalArraysToCalculate()
{
    // Check if the elements need distance matrix and/or radial distance vector
    for (auto &p : m_waveFunctionElements) {
        int need = p->getGlobalArrayNeed();
        if (need == 1) {
            m_calculateDistanceMatrix = true;
        } else if (need == 2) {
            m_calculateRadialVector = true;
        } else if (need == 3) {
            m_calculateDistanceMatrix = true;
            m_calculateRadialVector = true;
        }
    }
    // Check if the Hamiltonian needs distance matrix and/or radial distance vector
    int need = m_hamiltonian->getGlobalArrayNeed();
    if (need == 1) {
        m_calculateDistanceMatrix = true;
    } else if (need == 2) {
        m_calculateRadialVector = true;
    } else if (need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
    need = m_interactionStyle->getGlobalArrayNeed();
    if (need == 1) {
        m_calculateDistanceMatrix = true;
    } else if (need == 2) {
        m_calculateRadialVector = true;
    } else if (need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
}


/* ----------------------------------------------------------------------------
  Set the maximum number of parameters found in a wave function element. This
  will be used to create the parameter matrix, with dim
  (maxNumberOfParameter x numberOfElements). Determines also the total number
  of particles. Called automatically when needed
---------------------------------------------------------------------------- */

void System::setMaxParameters()
{
    int maxNumberOfElements = 0;
    int counter = 0;
    for (auto &i : m_waveFunctionElements) {
        int numberOfParameters = i->getNumberOfParameters();
        if (numberOfParameters > maxNumberOfElements) {
            maxNumberOfElements = numberOfParameters;
        }
        counter += numberOfParameters;
    }
    m_maxParameters = maxNumberOfElements;
    m_totalNumberOfParameters = counter;
}

