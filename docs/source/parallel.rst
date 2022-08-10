Parallelisation
================

The code supports parallelisation through MPI, where the system is sampled independently on each process. This requires just a few KB of communication before and after the parameters are updated, with a time neglectable compared to the sampling time. However, the burn-in time should remain the same no matter how many processes that is used. To run the code on 4 processes, execute 

.. code-block:: bash

   mpirun -n 4 vmachine input.in

The code has been running on thousands of processes in parallel without any problem. 
