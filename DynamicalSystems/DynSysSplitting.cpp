//
// Created by Nikita Kruk on 29.06.20.
//

#include "DynSysSplitting.hpp"

#include <mpi.h>
#include <algorithm> // std::min
#include <cassert>

DynSysSplitting::DynSysSplitting(Thread *thread) :
    thread_(thread),
    sin_of_phase_(kL, 0.0),
    cos_of_phase_(kL, 0.0),
//    density_slope_(kL * kK, 0.0),
    density_slope_wrt_angle_(kL * kK, 0.0),
    density_slope_wrt_angular_velocity_(kL * kK, 0.0),
    flux_limiter_count_(kDim, 0)
{
  for (int k = 0; k < kL; ++k)
  {
    sin_of_phase_[k] = std::sin(Phi(k));
    cos_of_phase_[k] = std::cos(Phi(k));
  } // k
}

DynSysSplitting::~DynSysSplitting()
{
  sin_of_phase_.clear();
  cos_of_phase_.clear();
//  density_slope_.clear();
  density_slope_wrt_angle_.clear();
  density_slope_wrt_angular_velocity_.clear();
}

Real DynSysSplitting::CalculateAngularUpdate(const std::vector<Real> &system_state,
                                             const std::vector<Real> &k_prev,
                                             std::vector<Real> &k_next,
                                             Real k_coef,
                                             Real dt)
{
  cfl_a_ = 0.0;

  static std::vector<Real> rk_system_state(system_state.size(), 0.0);
  for (int idx = 0, size = system_state.size(); idx < size; ++idx)
  {
    rk_system_state[idx] = system_state[idx] + k_coef * dt * k_prev[idx];
  } // idx

  static std::vector<Real> angular_flux(system_state.size(), (Real) 0.0);
  std::fill(angular_flux.begin(), angular_flux.end(), Real(0.0));
//  CalculateFlux(Dimension::kPhi, rk_system_state, angular_flux);
  CalculateAngularFlux(system_state, angular_flux);
  thread_->SynchronizeVector(angular_flux);

  int i = 0, j = 0;
  int idx_prev_phi = 0;
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev_phi);
    k_next[idx] = -(angular_flux[idx] - angular_flux[idx_prev_phi]) / kDphi;
  } // idx

  /*if (thread_->IsRoot())
  {
    std::cout << "CFL: " << dt << " <= " << std::min(kDx / (4.0 * cfl_a_), kDy / (4.0 * cfl_b_)) << std::endl;
  }*/
  return kDphi / (2.0 * cfl_a_);
}

Real DynSysSplitting::CalculateAngularVelocityUpdate(const std::vector<Real> &system_state,
                                                     const std::vector<Real> &k_prev,
                                                     std::vector<Real> &k_next,
                                                     Real k_coef,
                                                     Real dt)
{
  cfl_b_ = 0.0;

  static std::vector<Real> rk_system_state(system_state.size(), 0.0);
  for (int idx = 0, size = system_state.size(); idx < size; ++idx)
  {
    rk_system_state[idx] = system_state[idx] + k_coef * dt * k_prev[idx];
  } // idx

  static std::vector<Real> angular_velocity_flux(system_state.size(), (Real) 0.0);
  std::fill(angular_velocity_flux.begin(), angular_velocity_flux.end(), Real(0.0));
//  CalculateFlux(Dimension::kOmega, rk_system_state, angular_velocity_flux);
  CalculateAngularVelocityFlux(system_state, angular_velocity_flux);
  thread_->SynchronizeVector(angular_velocity_flux);

  int i = 0, j = 0;
  int idx_prev_omega = 0;
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev_omega);
    // impose zero flux boundary conditions
    if (j == 0)
    {
      k_next[idx] = -(angular_velocity_flux[idx] - 0.0) / kDomega;
    } else if (j == kK - 1)
    {
      k_next[idx] = -(0.0 - angular_velocity_flux[idx_prev_omega]) / kDomega;
    } else
    {
      k_next[idx] = -(angular_velocity_flux[idx] - angular_velocity_flux[idx_prev_omega]) / kDomega;
    }
  } // idx

  /*if (thread_->IsRoot())
  {
    std::cout << "CFL: " << dt << " <= " << kDphi / (2.0 * cfl_c_) << std::endl;
  }*/
  return kDphi / (2.0 * cfl_b_);
}

/*void DynSysSplitting::CalculateFlux(Dimension dim,
                                    const std::vector<Real> &system_state,
                                    std::vector<Real> &flux)
{

}*/

void DynSysSplitting::CalculateAngularFlux(const std::vector<Real> &system_state, std::vector<Real> &angular_flux)
{

}

void DynSysSplitting::CalculateAngularVelocityFlux(const std::vector<Real> &system_state,
                                                   std::vector<Real> &angular_velocity_flux)
{

}

/**
 *
 * @param dim : 0, 1, 2
 * @param system_state
 */
void DynSysSplitting::CalculateDensitySlopes(Dimension dim, const std::vector<Real> &system_state)
{
  int i = 0, j = 0;
  int idx_next = 0, idx_prev = 0;
  Real delta = 0.0;
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    switch (dim)
    {
      case Dimension::kPhi : utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
        utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev);
        delta = kDphi;
        density_slope_wrt_angle_[idx] = (system_state[idx_next] - system_state[idx_prev]) / (2.0 * delta);
        break;
      case Dimension::kOmega : utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
        utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev);
        delta = kDomega;
        density_slope_wrt_angular_velocity_[idx] = (system_state[idx_next] - system_state[idx_prev]) / (2.0 * delta);
        break;
    }
//    density_slope_[idx] = (system_state[idx_next] - system_state[idx_prev]) / (2.0 * delta);
  } // idx
}

/**
 *
 * @param dim : 0, 1, 2
 * @param system_state
 */
void DynSysSplitting::VerifyPositivityOfDensityAtCellInterfaces(Dimension dim,
                                                                const std::vector<Real> &system_state,
                                                                const std::vector<Real> &density_slope)
{
  Real density_next = 0.0, density_prev = 0.0;
  Real delta = 0.0;
  switch (dim)
  {
    case Dimension::kPhi : delta = kDphi;
      break;
    case Dimension::kOmega : delta = kDomega;
      break;
  }
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    density_next = system_state[idx] + delta / 2.0 * density_slope[idx];
    density_prev = system_state[idx] - delta / 2.0 * density_slope[idx];
//    if (density_prev < 0.0)
//    {
    if ((density_prev < 0.0) || (density_next < 0.0))
    {
      RecalculateDensitySlopeWithFluxLimiter(dim, system_state, idx);
    }
//    } else if (density_next < 0.0)
//    {
//      RecalculateDensitySlopeWithFluxLimiter(dim, system_state, idx);
//    }
  } // idx
}

void DynSysSplitting::RecalculateDensitySlopeWithFluxLimiter(Dimension dim,
                                                             const std::vector<Real> &system_state,
                                                             int idx)
{
  int i = 0, j = 0, k = 0;
  int idx_next = 0, idx_prev = 0;
  Real delta = 0.0;
  utilities::OneDimIdxToTwoDimIdx(idx, i, j);

  const Real theta = 1.0;
  switch (dim)
  {
    case Dimension::kPhi : utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
      utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev);
      delta = kDphi;
      density_slope_wrt_angle_[idx] =
          utilities::Minmod(theta * (system_state[idx_next] - system_state[idx]) / delta,
                            density_slope_wrt_angle_[idx],
                            theta * (system_state[idx] - system_state[idx_prev]) / delta);
      break;
    case Dimension::kOmega : utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
      utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev);
      delta = kDomega;
      density_slope_wrt_angular_velocity_[idx] =
          utilities::Minmod(theta * (system_state[idx_next] - system_state[idx]) / delta,
                            density_slope_wrt_angular_velocity_[idx],
                            theta * (system_state[idx] - system_state[idx_prev]) / delta);
      break;
  }

//  density_slope_[idx] = utilities::Minmod(theta * (system_state[idx_next] - system_state[idx]) / delta,
//                                          density_slope_[idx],
//                                          theta * (system_state[idx] - system_state[idx_prev]) / delta);
  /*// max possible density slope
  if (left_is_negative)
  {
    density_slope_[idx] = system_state[idx] / (delta / 2.0);
  } else
  {
    density_slope_[idx] = -system_state[idx] / (delta / 2.0);
  }*/

  switch (dim)
  {
    case Dimension::kPhi : ++flux_limiter_count_[0];
      break;
    case Dimension::kOmega : ++flux_limiter_count_[1];
      break;
  }
}

void DynSysSplitting::ResetFluxLimiterCount()
{
  std::fill(flux_limiter_count_.begin(), flux_limiter_count_.end(), 0);
}