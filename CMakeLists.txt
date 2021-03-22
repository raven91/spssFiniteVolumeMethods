cmake_minimum_required(VERSION 3.16)
project(spssFiniteVolumeMethods)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include/)
link_directories(/usr/local/lib/)

set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -O3 -flto")
set(GCC_COVERAGE_LINK_FLAGS "-lmpc -lmpfr -lgmp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

add_executable(spssFiniteVolumeMethods
        main.cpp
        Definitions.hpp Definitions.cpp
        Parallelization/Parallelization.hpp Parallelization/Parallelization.cpp
        Parallelization/Thread.cpp Parallelization/Thread.hpp
        Parallelization/ThreadOmegaSharedMemory.cpp Parallelization/ThreadOmegaSharedMemory.hpp
        ControlEngines/InitialConditions.hpp
        ControlEngines/SimulationEnginePtr.cpp ControlEngines/SimulationEnginePtr.hpp
        Observers/BinaryObserverPtr.cpp Observers/BinaryObserverPtr.hpp
        Steppers/RungeKutta2StepperWithSplittingPtr.cpp Steppers/RungeKutta2StepperWithSplittingPtr.hpp
        DynamicalSystems/DynSysSplittingPtr.cpp DynamicalSystems/DynSysSplittingPtr.hpp
        DynamicalSystems/DynSysSplittingFastPtr.cpp DynamicalSystems/DynSysSplittingFastPtr.hpp
        FourierTransforms/Convolution.cpp FourierTransforms/Convolution.hpp
        FourierTransforms/FastDiscreteInverseFourierTransform.cpp FourierTransforms/FastDiscreteInverseFourierTransform.hpp
        ControlEngines/SimulationEngine.cpp ControlEngines/SimulationEngine.hpp
        Observers/BinaryObserver.cpp Observers/BinaryObserver.hpp
        Steppers/RungeKutta2StepperWithSplitting.cpp Steppers/RungeKutta2StepperWithSplitting.hpp
        DynamicalSystems/DynSysSplitting.cpp DynamicalSystems/DynSysSplitting.hpp
        DynamicalSystems/DynSysSplittingFast.cpp DynamicalSystems/DynSysSplittingFast.hpp
        Parallelization/ThreadRandom.cpp Parallelization/ThreadRandom.hpp
        DynamicalSystems/DynSysSplittingRandomPtr.cpp DynamicalSystems/DynSysSplittingRandomPtr.hpp
        ControlEngines/SimulationEngineRandomPtr.cpp ControlEngines/SimulationEngineRandomPtr.hpp)