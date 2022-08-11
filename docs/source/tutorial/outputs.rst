Outputs
========

The code is able to write several different output files with various purposes. The log that is printed to the terminal contains all information about the simulation, and should be kept if one wants to redo the simulation. To write the log file, use

.. code-block::

   vmachine input.in >> vmachine.log

Energies
---------

Sometimes one wants to inspect the energy expectation values during the simulation, for instance to check if the energy seems to have converged. To write the energy expectation values to file, use

.. code-block::

    printEnergyToFile: true

Then both the total energy (energy.dat), kinetic energy (kinetic.dat), oscillator energy (external.dat) and interaction energy (interaction.dat) is written to file. 

Electron densities
--------------------

Electron densities is often a convenient way of representing the wave function and mapping the high-dimensional space down to lower dimensions. The code supports radial and spatial onebody density and radial twobody density. To write these out to file, use

.. code-block::

   computeOneBodyDensity: true
   computeOneBodyDensity2: true
   computeTwoBodyDensity: true

for the respective electron densities. Files containing the densities on the grid are then written to onebody.dat, onebody2.dat and twobody.dat. To specify the number of bins and the maximum range of the grid (maximum radius), use

.. code-block::

   numberOfBins: 100
   maxRadius: 5

Parameters
-----------

The parameters are usually written out to file for checkpointing purposes, but the parameters may also be analysed. To print parameters to file, use

.. code-block::

   printParametersToFile: true

The checkpoint frequency can be set using

.. code-block::

   checkpointFreq: 100

The default checkpoint frequency is 100.
