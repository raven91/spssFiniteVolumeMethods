//
// Created by Nikita Kruk on 29.06.20.
//

#include "SimulationEngineRandomPtr.hpp"
#include "InitialConditions.hpp"
#include "../Observers/BinaryObserverPtr.hpp"
#include "../Steppers/RungeKutta2StepperWithSplittingPtr.hpp"
#include "../DynamicalSystems/DynSysSplittingRandomPtr.hpp"

SimulationEngineRandomPtr::SimulationEngineRandomPtr(ThreadRandom *thread) :
    thread_(thread),
    system_state_size_(kL * kK)
{
  system_state_ = new Real[system_state_size_];
}

SimulationEngineRandomPtr::~SimulationEngineRandomPtr()
{
  delete[] system_state_;
}

void SimulationEngineRandomPtr::InitializeSystemStateByRule()
{
  if (thread_->IsRoot())
  {
    InitialConditions::UniformPerturbed < Real *const>(system_state_, system_state_size_);
  }

  // normalize the probability density function
  Real overall_density =
      std::accumulate(&system_state_[0], &system_state_[system_state_size_], Real(0.0)) * kDphi * kDomega;
  std::for_each(&system_state_[0], &system_state_[system_state_size_], [&](Real &ss) { ss /= overall_density; });

  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEngineRandomPtr::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  InitializeSystemStateByRule();
  BinaryObserverPtr observer(thread_);
  DynSysSplittingRandomPtr system(thread_);
  RungeKutta2StepperWithSplittingPtr stepper(thread_, kDt);

  Real t = kT0;
  observer.SaveSystemState(system_state_, system_state_size_, t);
  while (t <= kT1)
  {
    t += stepper.GetDt();
    stepper.DoStep(system, system_state_);
    observer.SaveSystemState(system_state_, system_state_size_, t);
    observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
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

  if (thread_->IsRoot())
  {
    std::cout << "simulation ended" << std::endl;
  }
}