/* ----------------------------------------------------------------------------
  Mersenne-Twister pseudo random number generator (prng). The prng has a 
  period of 2^19937, and is popularly abbreviated mt19937. The algorithm was
  first proposed and described by ... This code is just a wrapper around the
  standard C++ implementation of Mersenne-Twister.
---------------------------------------------------------------------------- */

#include <iostream>
#include <random>

#include "mersennetwister.h"


/* ----------------------------------------------------------------------------
   Mersenne Twister constructor
---------------------------------------------------------------------------- */

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator(), generator(seed())
{}

//Mersenne Twister RNG
//std::random_device rd;  //Will be used to obtain a seed for the random number engine
//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

/* ----------------------------------------------------------------------------
   Set seed
---------------------------------------------------------------------------- */

void MersenneTwister::setSeed(int seed)
{
    generator.seed(seed);
}


/* ----------------------------------------------------------------------------
   Returns a double drawn from a Gaussian distribution with mean 'mean' and
   variance 'variance'
---------------------------------------------------------------------------- */

double MersenneTwister::nextGaussian(double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return dis(generator);
}


/* ----------------------------------------------------------------------------
   Returns an interger between 0 and (including) 'upperLimit'
---------------------------------------------------------------------------- */

int MersenneTwister::nextInt(int upperLimit)
{
    std::uniform_int_distribution<int> dis(0, upperLimit - 1);
    return dis(generator);
}


/* ----------------------------------------------------------------------------
   Returns a double between 0 and 1 drawn from a uniform distribution
---------------------------------------------------------------------------- */

double MersenneTwister::nextDouble(double low, double high)
{
    std::uniform_real_distribution<double> dis(low, high);
    return dis(generator);
}

Eigen::MatrixXd MersenneTwister::randomUniformMatrix(Eigen::Index row, Eigen::Index col, double low, double high)
{
    std::uniform_real_distribution<double> dis(low, high);
    return Eigen::MatrixXd::Zero(row, col).unaryExpr([&](double /*dummy*/){return dis(generator);});
}

Eigen::MatrixXd MersenneTwister::randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return Eigen::MatrixXd::Zero(row, col).unaryExpr([&](double /*dummy*/){return dis(generator);});
}
