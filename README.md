# VMaChine
VMaChine is a general variational Monte Carlo (VMC) solver written in object-oriented C++. It was implemented with focus on quantum dot systems, and standard trial wave functions (Hermite functions) are implemented. However, we primary aim of examining trial wave functions where a less amount of physical intuition is required. Therefore, several trial wave functions based on neural networks and machine learning are also included, hence VMaChine. The code is also largely parallelizable and is made to be fast.

## Basic usage
After compiling the code, a configuration file ```input.in``` is run by
```bash
vmachine input.in
```
on a single thread or
```bash
mpirun -n 4 vmachine input.in
```
on 4 CPUs (might be any number).

## What's new?

#### Update 2021-01-01
- An interaction class is added. The old style "interaction: true" is now superseded by "interactionStyle: coulomb".

#### Update 2020-12-28
- The log (printed to the terminal by default) is now more readable both for a human eye and for the computer
- Examples added

#### Update 2020-12-27
- Fixed bug in CMakeLists.txt
- Removed support for Armadillo, as it was found to be more confusing than useful


## Prerequisites
To run the code without issues, C++14 is required. In addition, a few external packages are needed:
- MPI
- Eigen
- Blocker

#### MPI
MPI is used for parallel processing. On Linux, the package can be installed by the following commands
```bash
sudo apt-get install libopenmpi-dev
sudo apt-get install openmpi-bin
```
MPI is also avaliable on other platforms.

#### Eigen
Eigen is a C++ template library for linear algebra operations. See
[http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/) for installation details.

#### Blocker
Blocker is an auto blocking package developed by Marius Jonsson, which is our preferred resampling tool. To get the package, go to [https://github.com/computative/block](https://github.com/computative/block) and clone the repository:

```bash
git clone https://github.com/computative/block.git
```

## Install
VMaChine can be installed by cloning this repository. This is preferably done in the home directory:
```bash
cd
git clone https://github.com/evenmn/VMaChine.git
```
Then copy Eigen header files and blocker files to VMaChine:
```bash
cd ~/Download/eigen-3.3.9   # insert correct path here
cp -r Eigen ~/VMaChine/src/
mv ~/block ~/VMaChine/src/
```

## Build
The code can be compiled by either CMake or QMake

#### CMake (recommended)
Build VMaChine in the usual CMake way:
```bash
mkdir build
cd build
cmake ..
make -j8
```
The executable is then found in ```~/VMaChine/build``` folder, which is added to the bash resource path by
```bash
echo "export PATH=~/VMaChine/build:\$PATH" >> ~/.bashrc
```

#### QMake (QT-creator)
1. [Download QT-creator](https://www.qt.io/download-qt-installer?hsCtaTracking=9f6a2170-a938-42df-a8e2-a9f0b1d6cdce%7C6cb0de4f-9bb5-4778-ab02-bfb62735f3e5)
2. Configure the building file ```src/vmachine.pro```

The project can then be run in QT-creator using ```ctrl``` + ```R```.

-------------------

## Running simulation
The desired simulation parameters are specified in a configuration file. A typical configuration file can be found in [input.in](examples/quantumdot_vmc/input.in):

```bash
# system
numParticles: 2
numDimensions: 2
hamiltonian: harmonicOscillator
omega: 1.0
interactionStyle: coulomb

# wave function
waveFunctionElement: gaussian
waveFunctionElement: padeJastrow

# simulation
numIterations: 100
numSteps: 100000
learningRate: 0.1
stepLength: 0.1
```
The configuration file is simply run by
```bash
vmachine input.in
```

## Licence
[MIT](https://choosealicense.com/licenses/mit/)
