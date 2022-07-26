#pragma once
#include <iostream>
#include <random>
#include <Eigen/Dense>

class RandomNumberGenerator
{
public:
    RandomNumberGenerator();
    virtual void setSeed(int seed) = 0;
    virtual int nextInt(int upperLimit) = 0;
    virtual double nextDouble(double = 0.0, double = 1.0) = 0;
    virtual double nextGaussian(double mean, double variance) = 0;
    virtual std::string getLabel() = 0;
    virtual Eigen::MatrixXd randomUniformMatrix(Eigen::Index row, Eigen::Index col, double = 0.0, double = 1.0) = 0;
    virtual Eigen::MatrixXd randomNormalMatrix(Eigen::Index row, Eigen::Index col, double mean, double variance) = 0;

    virtual ~RandomNumberGenerator() = 0;

protected:
    std::random_device seed;
};
