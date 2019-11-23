#pragma once
#include <mpi.h>

struct parameters {
  double radius;
  double learningRate;
  double stepLength;
  int iterations;
};

parameters getParameters(int numberOfDimensions, int numberOfParticles, double omega) {

    // === D=2, N=2 ===

    if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 10.0) {
        parameters D2_N2_w10p0;
        D2_N2_w10p0.radius = 3.0;
        D2_N2_w10p0.learningRate = 0.5;
        D2_N2_w10p0.stepLength = 0.1;
        D2_N2_w10p0.iterations = 2000;
        return D2_N2_w10p0;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 5.0) {
        parameters D2_N2_w5p0;
        D2_N2_w5p0.radius = 3.0;
        D2_N2_w5p0.learningRate = 0.5;
        D2_N2_w5p0.stepLength = 0.1;
        D2_N2_w5p0.iterations = 2000;
        return D2_N2_w5p0;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 3.0) {
        parameters D2_N2_w3p0;
        D2_N2_w3p0.radius = 3.0;
        D2_N2_w3p0.learningRate = 0.5;
        D2_N2_w3p0.stepLength = 0.1;
        D2_N2_w3p0.iterations = 2000;
        return D2_N2_w3p0;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 2.0) {
        parameters D2_N2_w2p0;
        D2_N2_w2p0.radius = 3.0;
        D2_N2_w2p0.learningRate = 0.5;
        D2_N2_w2p0.stepLength = 0.1;
        D2_N2_w2p0.iterations = 2000;
        return D2_N2_w2p0;
    }

    else if (numberOfDimensions == 2 && numberOfParticles == 2 && omega == 1.0) {
        parameters D2_N2_w1p0;
        D2_N2_w1p0.radius = 3.0;
        D2_N2_w1p0.learningRate = 0.5;
        D2_N2_w1p0.stepLength = 0.1;
        D2_N2_w1p0.iterations = 2000;
        return D2_N2_w1p0;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 0.5) {
        parameters D2_N2_w0p5;
        D2_N2_w0p5.radius = 4.0;
        D2_N2_w0p5.learningRate = 0.5;
        D2_N2_w0p5.stepLength = 0.1;
        D2_N2_w0p5.iterations = 2000;
        return D2_N2_w0p5;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 0.28) {
        parameters D2_N2_w0p28;
        D2_N2_w0p28.radius = 5.0;
        D2_N2_w0p28.learningRate = 0.5;
        D2_N2_w0p28.stepLength = 0.1;
        D2_N2_w0p28.iterations = 2000;
        return D2_N2_w0p28;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 0.1) {
        parameters D2_N2_w0p28;
        D2_N2_w0p28.radius = 5.0;
        D2_N2_w0p28.learningRate = 0.5;
        D2_N2_w0p28.stepLength = 0.1;
        D2_N2_w0p28.iterations = 2000;
        return D2_N2_w0p28;
    }

    else if(numberOfDimensions== 2 && numberOfParticles == 2 && omega == 0.01) {
        parameters D2_N2_w0p28;
        D2_N2_w0p28.radius = 5.0;
        D2_N2_w0p28.learningRate = 0.5;
        D2_N2_w0p28.stepLength = 0.1;
        D2_N2_w0p28.iterations = 2000;
        return D2_N2_w0p28;
    }

    else {
        perror("Required hyper-parameters are not specified and getting the automatized parameters failed");
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
}
