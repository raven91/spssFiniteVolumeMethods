//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINERANDOMPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINERANDOMPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadRandom.hpp"

class SimulationEngineRandomPtr
{
 public:

  explicit SimulationEngineRandomPtr(ThreadRandom *thread);
  ~SimulationEngineRandomPtr();

  void RunSimulation();

 private:

  ThreadRandom *thread_;
  Real *system_state_;
  long system_state_size_;

  void InitializeSystemStateByRule();

};

#endif //SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINERANDOMPTR_HPP
