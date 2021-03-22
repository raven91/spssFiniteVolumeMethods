//
// Created by Nikita Kruk on 29.06.20.
//

#include "SimulationEngine.hpp"
#include "InitialConditions.hpp"
#include "../Observers/BinaryObserver.hpp"
#include "../Steppers/RungeKutta2StepperWithSplitting.hpp"
#include "../DynamicalSystems/DynSysSplittingFast.hpp"

SimulationEngine::SimulationEngine(Thread *thread) :
    thread_(thread),
    system_state_((unsigned long) (kL * kK), 0.0)
{

}

SimulationEngine::~SimulationEngine()
{
  system_state_.clear();
}

void SimulationEngine::InitializeSystemStateByRule()
{
  if (thread_->IsRoot())
  {
    InitialConditions::UniformPerturbed<std::vector<Real>>(system_state_, system_state_.size());

    // normalize the probability density function
    Real overall_density =
        std::accumulate(system_state_.begin(), system_state_.end(), Real(0.0)) * kDphi * kDomega;
    std::for_each(system_state_.begin(), system_state_.end(), [&](Real &ss) { ss /= overall_density; });
  }

  thread_->BroadcastVector(system_state_);
}

#if defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
#endif
void SimulationEngine::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

#if defined(MPI_PARAMETER_SCAN)
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  kAlpha = 0.0 + rank * 0.1;
#endif

  InitializeSystemStateByRule();
  BinaryObserver observer(thread_);
  DynSysSplittingFast system(thread_);
  RungeKutta2StepperWithSplitting stepper(thread_, kDt);

  Real t = kT0;
  observer.SaveSystemState(system_state_, t);
  while (t <= kT1)
  {
    t += stepper.GetDt();
    stepper.DoStep(system, system_state_);
    observer.SaveSystemState(system_state_, t);
    observer.SaveSummaryStatistics(system_state_, t);
    observer.SaveAdditionalInformation(stepper.GetAverageFluxLimiterCount(),
                                       system.GetCflA(),
                                       system.GetCflB(),
                                       t);
    thread_->BroadcastCondition(observer.should_terminate_);
    if (observer.should_terminate_)
    {
      break;
    }
  }
#if defined(MPI_PARAMETER_SCAN)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (thread_->IsRoot())
  {
    std::cout << "simulation ended" << std::endl;
  }
}