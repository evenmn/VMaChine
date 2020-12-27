#include "mersennetwister.h"
#include "../system.h"

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator()
{}

//Mersenne Twister RNG
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

double MersenneTwister::nextGaussian(double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return dis(gen);
}

int MersenneTwister::nextInt(int upperLimit)
{
    std::uniform_int_distribution<int> dis(0, upperLimit - 1);
    return dis(gen);
}

double MersenneTwister::nextDouble()
{
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}

Eigen::MatrixXd MersenneTwister::randomUniformMatrix(Eigen::Index row, Eigen::Index col)
{
    std::uniform_real_distribution<double> dis(0, 1);
    return Eigen::MatrixXd::Zero(row, col).unaryExpr([&](double /*dummy*/){return dis(gen);});
}

Eigen::MatrixXd MersenneTwister::randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance)
{
    std::normal_distribution<double> dis(mean, variance);
    return Eigen::MatrixXd::Zero(row, col).unaryExpr([&](double /*dummy*/){return dis(gen);});
}
