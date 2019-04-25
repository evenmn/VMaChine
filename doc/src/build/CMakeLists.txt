project(QMC)

# Add directories containing header-files
include_directories(
    .
    Basis
    Hamiltonians
    InitialStates
    InitialWeights
    Metropolis
    Optimization
    Resampling
    RNG
    WaveFunctions
)

# Add sources from all directories
file(
    GLOB SOURCES "*.cpp"
    GLOB SOURCES "Basis/*.cpp"
    GLOB SOURCES "Hamiltonians/*.cpp"
    GLOB SOURCES "InitialStates/*.cpp"
    GLOB SOURCES "InitialWeights/*.cpp"
    GLOB SOURCES "Metropolis/*.cpp"
    GLOB SOURCES "Optimization/*.cpp"
    GLOB SOURCES "Resampling/*.cpp"
    GLOB SOURCES "RNG/*.cpp"
    GLOB SOURCES "WaveFunctions/*.cpp"
)


# Create executable called "vmc"
add_executable(vmc ${SOURCES})
