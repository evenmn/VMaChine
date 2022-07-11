#pragma once
#include <random>
#include <Eigen/Dense>

#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    void setSeed(int seed) override;
    int nextInt(int upperLimit) override;
    double nextDouble() override;
    double nextGaussian(double mean, double variance) override;
    std::string getLabel() override { return m_label; }
    Eigen::MatrixXd randomUniformMatrix(Eigen::Index row, Eigen::Index col);
    Eigen::MatrixXd randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance);
private:
    std::mt19937 generator;
    std::string m_label = "Mersenne-Twister";
};
