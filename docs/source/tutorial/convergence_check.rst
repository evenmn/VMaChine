Convergence Check
==================

The implemented convergence check is rather simple: If the energy does not vary with more than some tolerance over a given range of energies, the simulation has converged. This is checked using

.. code-block::

   checkConvergence: true

The tolerance is set by:

.. code-block::

   tolerance: 0.001

And the energy range is given by

.. code-block::

   numberOfEnergies: 10


