//
// Created by Nikita Kruk on 29.06.20.
//

#include "RungeKutta2StepperWithSplitting.hpp"

#include <algorithm> // std::min_element

RungeKutta2StepperWithSplitting::RungeKutta2StepperWithSplitting(Thread *thread, Real dt) :
    thread_(thread),
    dt_(dt),
    order_(2),
    k_1_((unsigned long) (kL * kK), 0.0),
    k_2_((unsigned long) (kL * kK), 0.0),
    average_flux_limiter_count_(kDim, 0)
{

}

RungeKutta2StepperWithSplitting::~RungeKutta2StepperWithSplitting()
{
  k_1_.clear();
  k_2_.clear();
}

void RungeKutta2StepperWithSplitting::DoStep(DynSysSplitting &system, std::vector<Real> &system_state)
{
// Three-step second order method
  system.ResetFluxLimiterCount();

  Real new_dts[3] = {0.0};
  new_dts[0] = DoAngularStep(system, system_state, dt_ / 2.0);
  new_dts[0] *= 2.0;

  new_dts[1] = DoAngularVelocityStep(system, system_state, dt_);

  new_dts[2] = DoAngularStep(system, system_state, dt_ / 2.0);
  new_dts[2] *= 2.0;

  FindAveragedFluxLimiterCounts(system);
//  dt_ = *std::min_element(new_dts, new_dts + 3);
//  thread_->FindMinValue(dt_);
}

Real RungeKutta2StepperWithSplitting::DoAngularStep(DynSysSplitting &system, std::vector<Real> &system_state, Real dt)
{
  std::vector<Real> new_dts(order_, 0.0);

  // No need to reset k_1_ to 0.0, since k_coef is 0.0
  new_dts[0] = system.CalculateAngularUpdate(system_state, k_1_, k_1_, 0.0, dt);
  thread_->SynchronizeVector(k_1_);
  new_dts[1] = system.CalculateAngularUpdate(system_state, k_1_, k_2_, 1.0, dt);

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    system_state[idx] += (k_1_[idx] + k_2_[idx]) * dt / 2.0;
  } // idx
  thread_->SynchronizeVector(system_state);

  return *std::min_element(new_dts.begin(), new_dts.end());
}

Real RungeKutta2StepperWithSplitting::DoAngularVelocityStep(DynSysSplitting &system,
                                                            std::vector<Real> &system_state,
                                                               Real dt)
{
  std::vector<Real> new_dts(order_, 0.0);

  // No need to reset k_1_ to 0.0, since k_coef is 0.0
  new_dts[0] = system.CalculateAngularVelocityUpdate(system_state, k_1_, k_1_, 0.0, dt);
  thread_->SynchronizeVector(k_1_);
  new_dts[1] = system.CalculateAngularVelocityUpdate(system_state, k_1_, k_2_, 1.0, dt);

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    system_state[idx] += (k_1_[idx] + k_2_[idx]) * dt / 2.0;
  } // idx
  thread_->SynchronizeVector(system_state);

  return *std::min_element(new_dts.begin(), new_dts.end());
}

void RungeKutta2StepperWithSplitting::FindAveragedFluxLimiterCounts(const DynSysSplitting &system)
{
  average_flux_limiter_count_ = system.GetFluxLimiterCount();
//  if (thread_->IsRoot())
//  {
//    std::cout << "counts: " << average_flux_limiter_count_[0] << " " << average_flux_limiter_count_[1] << " "
//              << std::endl;
//  }
  thread_->ComputeSumIntoRootOnly(average_flux_limiter_count_);
  average_flux_limiter_count_[0] = average_flux_limiter_count_[0] / (2 * order_); // 2 for two advection steps
  average_flux_limiter_count_[1] = average_flux_limiter_count_[1] / (order_);
}