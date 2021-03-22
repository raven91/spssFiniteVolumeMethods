//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTINGPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTINGPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/DynSysSplittingPtr.hpp"
#include "../Parallelization/ThreadOmegaSharedMemory.hpp"

#include <vector>

class RungeKutta2StepperWithSplittingPtr
{
 public:

  explicit RungeKutta2StepperWithSplittingPtr(ThreadOmegaSharedMemory *thread, Real dt);
  ~RungeKutta2StepperWithSplittingPtr();

  [[nodiscard]] Real GetDt() const
  { return dt_; }
  [[nodiscard]] const std::vector<int> &GetAverageFluxLimiterCount() const
  { return average_flux_limiter_count_; }

  void DoStep(DynSysSplittingPtr &system, Real *const system_state);

 private:

  ThreadOmegaSharedMemory *thread_;
  Real dt_;
  int order_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<int> average_flux_limiter_count_;

  Real DoAngularStep(DynSysSplittingPtr &system,
                       Real *const system_state,
                       Real dt);
  Real DoAngularVelocityStep(DynSysSplittingPtr &system,
                             Real *const system_state,
                             Real dt);
  void FindAveragedFluxLimiterCounts(const DynSysSplittingPtr &system);

};

#endif //SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTINGPTR_HPP
