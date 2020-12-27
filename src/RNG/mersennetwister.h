#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int nextInt(int upperLimit);
    double nextDouble();
    double nextGaussian(double mean, double variance);
    Eigen::MatrixXd randomUniformMatrix(Eigen::Index row, Eigen::Index col);
    Eigen::MatrixXd randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance);
};
