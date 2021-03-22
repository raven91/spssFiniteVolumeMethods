//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTING_HPP
#define SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTING_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/DynSysSplitting.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class RungeKutta2StepperWithSplitting
{
 public:

  explicit RungeKutta2StepperWithSplitting(Thread *thread, Real dt);
  ~RungeKutta2StepperWithSplitting();

  [[nodiscard]] Real GetDt() const
  { return dt_; }
  [[nodiscard]] const std::vector<int> &GetAverageFluxLimiterCount() const
  { return average_flux_limiter_count_; }

  void DoStep(DynSysSplitting &system, std::vector<Real> &system_state);

 private:

  Thread *thread_;
  Real dt_;
  int order_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<int> average_flux_limiter_count_;

  Real DoAngularStep(DynSysSplitting &system, std::vector<Real> &system_state, Real dt);
  Real DoAngularVelocityStep(DynSysSplitting &system, std::vector<Real> &system_state, Real dt);
  void FindAveragedFluxLimiterCounts(const DynSysSplitting &system);

};

#endif //SPSSFINITEVOLUMEMETHODS_RUNGEKUTTA2STEPPERWITHSPLITTING_HPP
