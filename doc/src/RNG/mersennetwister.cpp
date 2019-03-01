#include "mersennetwister.h"
#include "rng.h"
#include <random>
#include <iostream>
#include "../system.h"

MersenneTwister::MersenneTwister()  :
    RandomNumberGenerator() {
}

//Mersenne Twister RNG
std::random_device rd;                       //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd());                      //Standard mersenne_twister_engine seeded with rd()


double MersenneTwister::nextGaussian(double mean, double standardDeviation) {
    std::normal_distribution<double> rand_gauss(mean, standardDeviation);
    return rand_gauss(gen);
}

int MersenneTwister::nextInt(int upperLimit) {
    std::uniform_int_distribution<> rand_int(0, upperLimit - 1);
    return rand_int(gen);
}

double MersenneTwister::nextDouble() {
    std::uniform_real_distribution<> rand_double(0, 1);
    return rand_double(gen);
}
