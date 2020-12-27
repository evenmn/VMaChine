# VMaChine
----------------------
VMaChine is a general variational Monte Carlo (VMC) solver written in object-oriented C++. It was implemented with focus on quantum dot systems, and standard trial wave functions (Hermite functions) are implemented. However, we primary aim of examining trial wave functions where a less amount of physical intuition is required. Therefore, several trial wave functions based on neural networks and machine learning are also included, hence VMaChine. The code is also largely parallelizable and is made to be fast.

## What's new?

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
[http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/) for installation details. Copy source directory to VMaChine:

```bash
cp -r Eigen ~/VMaChine/src/
```

#### Blocker
Blocker is an auto blocking package developed by Marius Jonsson, which is our preferred resampling tool. To get the package, go to [https://github.com/computative/block](https://github.com/computative/block) and clone the repository:

```bash
git clone https://github.com/computative/block.git
mv block ~/VMaChine/src/
```
-------------------

## Install
VMaChine can be installed by cloning this repository. This is preferably done in the home directory:
```bash
cd
git clone https://github.com/evenmn/VMaChine.git
```

-------------------

## Build
The code can be compiled by either CMake or QMake

#### CMake (recommended)
Build VMaChine in the standard CMake way:
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

##### Parallel processing using QT-creator
To run in parallel,  add an MPI executable. Go to ```Projects-> Run-> Add-> Custom Executable```. Then set
- Executable: ```/usr/bin/mpirun``` (or whatever ```which mpirun``` returns)
- Command line arguments: ```-n 4 VMaChine```
- Working directory: ```/where/the/executable/is```

This setup will run 4 parallel processes. The configuration window should look similar to this
![Run settings](screenshots/run_settings.png)

-------------------

## Setting up system
The parameters can be set in two different ways:

1. They can be specified in a configuration file (standard)
2. They can be set by editing ```main.cpp``` (requires the code to be recompiled every time changes are made, and is therefore usually not recommended)

In general, the settings from the config file overwrite the setting in ```main.cpp```, which again overwrites the default settings found in ```system.h```.

#### From configuration file
The configuration file, ```config```, needs to be passed as an argument when running the code,
```bash
mpirun -n 4 VMaChine config
```
The advantage of this is that we can change parameters without compiling the code again, which allows us to run multiple simulations in parallel. In particular, this is beneficial when running large simulations on computer clusters. To tell the program to access the config file, simply add
```c++
QD->initializeFromConfig(argc, argv);
```
to ```main.cpp```. An example on such a config file is found in [config](config).


#### From ```main.cpp```
From ```main.cpp```, all possible settings can be specified. Below, we present a minimalistic example on how it may look like for a two-dimensional quantum dot system with 6 electrons and frequency 1. We choose a Slater-Jastrow wave function as our trial wave function.
``` c++
#include "main.h"

int main(int argc, char** argv)
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(6);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setInteraction(true);

    QD->setLearningRate(0.1);
    QD->setStepLength(0.05);
    QD->setNumberOfMetropolisSteps(int(pow(2, 20)));

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function
    QD->setBasis(new Hermite(QD));
    QD->setWaveFunctionElement(new Gaussian(QD));
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    QD->setWaveFunctionElement(new PadeJastrow(QD));

    QD->setNumberOfIterations(1000);
    QD->runSimulation();
    return 0;
}
```
We first define the physical system and then we define method-related tools. ```main.h``` needs to be included in order to pass the objects ```HarmonicOscillator(QD)```, ```Hermite(QD)``` and so on. Apart from the system declaration and simulation call, the order of calls is irrelevant. Everything that is not specified in ```main.cpp``` will stay by the default.

-------------------

## Analyze results
### Energy and blocking results
The current energy is printed to the terminal for every iteration, together with the estimated variance, standard deviation, acceptence ratio and CPU time. _System info_ presents the settings used for the current run.

For the last iteration, blocking is performed and the blocking results are printed to the terminal in the very end. The terminal will typically look like this when a run is done:
![terminal](screenshots/screenshot_terminal.png)

The energies for all iterations are stored in a file found in the data folder. The exact location and file name is also printed to the terminal, see image above. To plot this energy file, run
```bash
python3 scripts/plot_energy.py
```
which also support multiple files.

### One-body density
The one-body density is calculated during the last iteration by default, and the file is stored in the same way as the energy files described above. The one-body density can be plotted by the command
```bash
python3 scripts/plot_ob_density.py
```
which again supports multiple files. Remember to set the correct number of dimensions inside the script!

-------------------

## Licence
[MIT](https://choosealicense.com/licenses/mit/)
