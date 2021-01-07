#pragma once
#include <iostream>
#include <random>
#include <Eigen/Dense>

class RandomNumberGenerator
{
public:
    RandomNumberGenerator();
    virtual int nextInt(int upperLimit) = 0;
    virtual double nextDouble() = 0;
    virtual double nextGaussian(double mean, double variance) = 0;
    virtual std::string getLabel() = 0;
    virtual Eigen::MatrixXd randomUniformMatrix(Eigen::Index row, Eigen::Index col) = 0;
    virtual Eigen::MatrixXd randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance) = 0;

    virtual ~RandomNumberGenerator() = 0;
};
