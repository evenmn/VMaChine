# VMC - General variational Monte-Carlo solver written in C++

## Requirements
To run this project without issues, the most recent C++ version, C++17, is recommended. In addition, a few external packages are required, that is
- Eigen (for linear algera operations)
- MPI (for parallel processing)

## Build & run
There are several ways to compile the code. Below, we will present two easy and rebust methods based on CMake and QMake. 

### CMake
1. chmod +x CompileVMC      # Giving CompileVMC the right permissions
2. ./CompileVMC             # Compile project
3. mpiexec -n 4 ./vmc       # Run code on 4 parallel processes
4. make clean               # Remove files 

### QMake (QT-creator)
1. [Download QT-creator](https://www.qt.io/download-qt-installer?hsCtaTracking=9f6a2170-a938-42df-a8e2-a9f0b1d6cdce%7C6cb0de4f-9bb5-4778-ab02-bfb62735f3e5)
2. Configurate the building file QMC.pro
3. Run using '''ctrl + R'''


## Adjust parameters
All necessary adjustments can currently be done in '''main.cpp'''
