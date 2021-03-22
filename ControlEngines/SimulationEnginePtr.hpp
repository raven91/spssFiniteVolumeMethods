//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINEPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINEPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadOmegaSharedMemory.hpp"

class SimulationEnginePtr
{
 public:

  explicit SimulationEnginePtr(ThreadOmegaSharedMemory *thread);
  ~SimulationEnginePtr();

  void RunSimulation();
  void RunContinuationMethod();

 private:

  ThreadOmegaSharedMemory *thread_;
  Real *system_state_;
  long system_state_size_;

  void InitializeSystemStateByRule();
  void InitializeSystemStateFromFile(bool should_add_perturbation);
  void AddPerturbationToSystemState();

};

#endif //SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINEPTR_HPP
