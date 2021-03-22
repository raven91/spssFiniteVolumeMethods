//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINE_HPP
#define SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class SimulationEngine
{
 public:

  explicit SimulationEngine(Thread *thread);
  ~SimulationEngine();

  void RunSimulation();

 private:

  Thread *thread_;
  std::vector<Real> system_state_;

  void InitializeSystemStateByRule();

};

#endif //SPSSFINITEVOLUMEMETHODS_SIMULATIONENGINE_HPP
