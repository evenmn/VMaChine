Installation
-------------

The code relies on two external C++ libraries, namely MPI and Eigen. It has been properly tested with OpenMPI only, but will probably also work with MPICH. On Ubuntu, install OpenMPI using

.. code-block:: bash

   sudo apt-get install libopenmpi-dev openmpi-bin

Download Eigen from `http://eigen.tuxfamily.org/ <http://eigen.tuxfamily.org/>`_ and put the :code:`Eigen` directory in :code:`VMaChine/include`.

Build
------

The code is based on the 2014 C++ standard, requiring g++ 6.1 or later. We recommend building the code with CMake. Follow the standard CMake instructions from top repo directory:

.. code-block:: bash

   mkdir build
   cd build
   cmake ..
   make -j8

The executable should appear in the :code:`build` directory. If something fails during the building process, please check g++ version and dependency links.
