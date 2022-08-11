Resampling
===========

The automatized blocking algorithm is implemented for resampling and estimate the variance of the energies. Both the algorithm itself and the code is written by Marius Jonsson. To perform resampling, simply set

.. code-block::

   doResampling: true

The blocking algorithm does only support data that is an integer power of 2 long, so make sure the number of steps is an integer power of 2. All sampling energies (both total, kinetic, external and interaction) will be written out to files with names starting with a random integer, which are (hopefully) deleted after the simulation. If the simulation is using several processes, they will write to separate files which is merged after sampling. 
