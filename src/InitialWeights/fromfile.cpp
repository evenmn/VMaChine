#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <cassert>

#include "fromfile.h"
#include "../system.h"

#define MAXBUFSIZE  ((int) 1e6)


/* ----------------------------------------------------------------------------
  Initialize weights from file class constructor
---------------------------------------------------------------------------- */

FromFile::FromFile(System *system, std::string filename)
    : InitialWeights(system)
{
    m_system = system;
    m_filename = filename;
}


/* ----------------------------------------------------------------------------
  Setup initial weights after number of elements and max number of parameters
  are found
---------------------------------------------------------------------------- */

void FromFile::setupInitialWeights()
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_parameters = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(m_filename);
    while (! infile.eof())
        {
        std::string line;
        std::getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    assert (rows == m_numberOfElements);
    assert (cols == m_maxParameters);

    // Populate matrix with numbers.
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m_parameters(i,j) = buff[ cols*i+j ];

    m_system->updateAllParameters(m_parameters);
}


/* ----------------------------------------------------------------------------
  Get parameters from class. Parameters are declared private
---------------------------------------------------------------------------- */

Eigen::MatrixXd FromFile::getParameters()
{
    return m_parameters;
}
