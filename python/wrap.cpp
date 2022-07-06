#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>        // automatic conversion between vector/valarray and numpy array
#include <pybind11/functional.h> // automatic conversion between functional and python function

#include "../src/RNG/mersennetwister.h"


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;


/* ----------------------------------------------------------------------------
   Python wrapper based on pybind11. Please see the pybind11 documentation
   for all function calls starting with "py::".
------------------------------------------------------------------------------- */

PYBIND11_MODULE(vmc, m) {

    // Mersenne Twister RNG
    py::class_<RandomNumberGenerator>(m, "RNG")
        .def("__repr__", &RandomNumberGenerator::getLabel,
            "Represent random number generator"
        );

    py::class_<MersenneTwister, RandomNumberGenerator>(m, "MersenneTwister")
        .def(py::init<>(),
            "Mersenne-Twister constructor"
        )
        .def("set_seed", &MersenneTwister::setSeed,
            "Set seed",
            py::arg("seed")
        )
        .def("next_int", &MersenneTwister::nextInt,
            "Get next integer in the sequence",
            py::arg("upper")
        )
        .def("next_double", &MersenneTwister::nextDouble,
            "Get next double (between 0 and 1) in the sequence"
        )
        .def("next_gaussian", &MersenneTwister::nextGaussian,
            "Get next Gaussian in the sequence",
            py::arg("mean") = 0.,
            py::arg("var") = 1.
        );
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
