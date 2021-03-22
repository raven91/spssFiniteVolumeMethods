#include "Parallelization/Parallelization.hpp"
#include "Parallelization/ThreadOmegaSharedMemory.hpp"
#include "Parallelization/ThreadRandom.hpp"
#include "ControlEngines/SimulationEngine.hpp"
#include "ControlEngines/SimulationEnginePtr.hpp"
#include "ControlEngines/SimulationEngineRandomPtr.hpp"

// Real kAlpha = 0.0;
//Real kDiffusionConstant = 0.0;

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);
  {
#if defined(MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY)
    ThreadOmegaSharedMemory thread(argc, argv);
    SimulationEnginePtr engine(&thread);
#elif defined(MPI_PARAMETER_SCAN)
    Thread thread(argc, argv);
    SimulationEngine engine(&thread);
#elif defined(MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM)
    ThreadRandom thread(argc, argv);
    SimulationEngineRandomPtr engine(&thread);
#else
    Thread thread(argc, argv);
    SimulationEngine engine(&thread);
#endif
    engine.RunSimulation();
//    engine.RunContinuationMethod();
  }
  FinalizeParallelSession();
  return 0;
}
