//
// Created by Nikita Kruk on 29.06.20.
//

#include "DynSysSplittingFast.hpp"
#include "../FourierTransforms/Convolution.hpp"

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <numeric> // std::accumulate
#include <algorithm> // std::fill

DynSysSplittingFast::DynSysSplittingFast(Thread *thread) :
    DynSysSplitting(thread),
    cos_kernel_at_spatial_point_(kL, 0.0),
    sin_kernel_at_spatial_point_(kL, 0.0)
{
  for (int i = 0; i < kL; ++i)
  {
    cos_kernel_at_spatial_point_[i] = (-kC2 + kC1) * std::cos(Phi(i) + kAlpha);
    sin_kernel_at_spatial_point_[i] = -kC1 * std::sin(Phi(i) + kAlpha);
  } // k
}

DynSysSplittingFast::~DynSysSplittingFast()
{

}

/*void DynSysSplittingFast::CalculateFlux(Dimension dim, const std::vector<Real> &system_state, std::vector<Real> &flux)
{
  Real density_next_in_current_cell = 0.0, density_prev_in_next_cell = 0.0;
  int idx = 0, idx_next = 0;
  Real delta = 0.0;
  switch (dim)
  {
    case Dimension::kPhi: delta = kDphi;
      break;
    case Dimension::kOmega: delta = kDomega;
      break;
  }

  CalculateDensitySlopes(dim, system_state);
  VerifyPositivityOfDensityAtCellInterfaces(dim, system_state);
  thread_->SynchronizeVector(density_slope_);

  static std::vector<Real> convolution(kL, 0.0);
  Real normalization = 0.0;
  CalculateConvolution(system_state, convolution, normalization);
  for (const int &angular_velocity_index: thread_->GetAngularVelocityLoopIndices())
  {
    int j = angular_velocity_index;
    for (int i = 0; i < kL; ++i)
    {
      utilities::TwoDimIdxToOneDimIdx(i, j, idx);
      density_next_in_current_cell = system_state[idx] + delta / 2.0 * density_slope_[idx];
      switch (dim)
      {
        case Dimension::kPhi: utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
          break;
        case Dimension::kOmega: utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
          break;
      }
      density_prev_in_next_cell = system_state[idx_next] - delta / 2.0 * density_slope_[idx_next];

      Real velocity = 0.0;
      switch (dim)
      {
        case Dimension::kPhi: velocity = Omega(j);
          if (cfl_a_ < std::fabs(velocity))
          {
            cfl_a_ = std::fabs(velocity);
          }
          break;
        case Dimension::kOmega:
          CalculateVelocityAtCellInterface(system_state,
                                           i,
                                           j,
                                           velocity,
                                           convolution,
                                           normalization);
          if (cfl_b_ < std::fabs(velocity))
          {
            cfl_b_ = std::fabs(velocity);
          }
          break;
      }
      flux[idx] = std::max((Real) 0.0, velocity) * density_next_in_current_cell
          + std::min((Real) 0.0, velocity) * density_prev_in_next_cell;
    } // i
  } // angular_velocity_index
}*/

void DynSysSplittingFast::CalculateAngularFlux(const std::vector<Real> &system_state, std::vector<Real> &angular_flux)
{
  Real density_next_in_current_cell = 0.0, density_prev_in_next_cell = 0.0;
  int idx_next = 0, i = 0, j = 0;
  Real velocity = 0.0;

  CalculateDensitySlopes(Dimension::kPhi, system_state);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kPhi, system_state, density_slope_wrt_angle_);
  thread_->SynchronizeVector(density_slope_wrt_angle_);

  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    velocity = Omega(j);
    if (velocity > 0.0)
    {
      density_next_in_current_cell = system_state[idx] + kDphi / 2.0 * density_slope_wrt_angle_[idx];
      angular_flux[idx] = velocity * density_next_in_current_cell;
    } else
    {
      utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
      density_prev_in_next_cell = system_state[idx_next] - kDphi / 2.0 * density_slope_wrt_angle_[idx_next];
      angular_flux[idx] = velocity * density_prev_in_next_cell;
    }
    if (cfl_a_ < std::fabs(velocity))
    {
      cfl_a_ = std::fabs(velocity);
    }
//    angular_flux[idx] = std::max((Real) 0.0, velocity) * density_next_in_current_cell
//        + std::min((Real) 0.0, velocity) * density_prev_in_next_cell;
//    if (velocity > 0.0)
//    {
//      angular_flux[idx] = velocity * density_next_in_current_cell;
//    } else
//    {
//      angular_flux[idx] = velocity * density_prev_in_next_cell;
//    }
  } // idx
}

void DynSysSplittingFast::CalculateAngularVelocityFlux(const std::vector<Real> &system_state,
                                                       std::vector<Real> &angular_velocity_flux)
{
  Real density_next_in_current_cell = 0.0, density_prev_in_next_cell = 0.0;
  int idx_next = 0, i = 0, j = 0;
  Real velocity = 0.0;

  // we need both density slopes for this flux
  CalculateDensitySlopes(Dimension::kPhi, system_state);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kPhi, system_state, density_slope_wrt_angle_);
  thread_->SynchronizeVector(density_slope_wrt_angle_);
  CalculateDensitySlopes(Dimension::kOmega, system_state);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kOmega, system_state, density_slope_wrt_angular_velocity_);
  thread_->SynchronizeVector(density_slope_wrt_angular_velocity_);
  // todo: synchronize only neighboring cell_omega

  static std::vector<Real> convolution(kL, 0.0);
  Real normalization = 0.0;
  CalculateConvolution(system_state, convolution, normalization);
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    CalculateVelocityAtCellInterface(system_state, i, j, velocity, convolution, normalization);
    if (velocity > 0.0)
    {
      density_next_in_current_cell = system_state[idx] + kDomega / 2.0 * density_slope_wrt_angular_velocity_[idx];
      angular_velocity_flux[idx] = velocity * density_next_in_current_cell;
    } else
    {
      utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
      density_prev_in_next_cell = system_state[idx_next] - kDomega / 2.0 * density_slope_wrt_angular_velocity_[idx_next];
      angular_velocity_flux[idx] = velocity * density_prev_in_next_cell;
    }
    if (cfl_b_ < std::fabs(velocity))
    {
      cfl_b_ = std::fabs(velocity);
    }
//    angular_velocity_flux[idx] = std::max((Real) 0.0, velocity) * density_next_in_current_cell
//        + std::min((Real) 0.0, velocity) * density_prev_in_next_cell;
//    if (velocity > 0.0)
//    {
//      angular_velocity_flux[idx] = velocity * density_next_in_current_cell;
//    } else
//    {
//      angular_velocity_flux[idx] = velocity * density_prev_in_next_cell;
//    }
  } // idx
}

/**
 * For the terms with convolution only
 * @param system_state
 * @param i
 * @param j
 * @param velocity
 */
void DynSysSplittingFast::CalculateVelocityAtCellInterface(const std::vector<Real> &system_state,
                                                           int i,
                                                           int j,
                                                           Real &velocity,
                                                           const std::vector<Real> &convolution,
                                                           Real normalization)
{
  Real potential_cur = 0.0, potential_next = 0.0;
  CalculateVelocityPotentialAtCellCenter(system_state, i, j, potential_cur, convolution, normalization);
  CalculateVelocityPotentialAtCellCenter(system_state,
                                         i,
                                         utilities::PositiveModulo(j + 1, kK),
                                         potential_next,
                                         convolution,
                                         normalization);
  velocity = -(potential_next - potential_cur) / kDomega;
}

/**
 * Calculation of a potential through convolution only
 * @param system_state
 * @param i
 * @param j
 * @param potentials
 */
void DynSysSplittingFast::CalculateVelocityPotentialAtCellCenter(const std::vector<Real> &system_state,
                                                                 int i,
                                                                 int j,
                                                                 Real &potential,
                                                                 const std::vector<Real> &convolution,
                                                                 Real normalization)
{
  Real omega_j = Omega(j);
  potential = 0.5 * kFriction * omega_j * omega_j - kSigma * omega_j * convolution[i] / normalization;
  if (kDiffusionConstant > Real(0))
  {
    int idx = 0;
    utilities::TwoDimIdxToOneDimIdx(i, j, idx);
    potential += (kDiffusionConstant * std::log(system_state[idx]));
  }
}

void DynSysSplittingFast::CalculateConvolution(const std::vector<Real> &system_state,
                                               std::vector<Real> &convolution,
                                               Real &normalization)
{
  static std::vector<Real> projected_system_state(kL, 0.0), projected_density_slope_wrt_angle(kL, 0.0);
  std::fill(projected_system_state.begin(), projected_system_state.end(), 0.0);
  std::fill(projected_density_slope_wrt_angle.begin(), projected_density_slope_wrt_angle.end(), 0.0);

  for (int m = 0; m < kK; ++m)
  {
    for (int l = 0; l < kL; ++l)
    {
      int idx = 0;
      utilities::TwoDimIdxToOneDimIdx(l, m, idx);
      projected_system_state[l] += system_state[idx];
      projected_density_slope_wrt_angle[l] += density_slope_wrt_angle_[idx];
    } // l
  } // m

  // convolution of the interaction kernel
  Convolution convolution_engine;
  static std::vector<Real> convolution_at_angular_velocity_point_0(kL, 0.0);
  convolution_engine.FastConvolution(projected_system_state,
                                     sin_kernel_at_spatial_point_,
                                     convolution_at_angular_velocity_point_0);
  static std::vector<Real> convolution_at_angular_velocity_point_1(kL, 0.0);
  convolution_engine.FastConvolution(projected_density_slope_wrt_angle,
                                     cos_kernel_at_spatial_point_,
                                     convolution_at_angular_velocity_point_1);

  for (int i = 0; i < kL; ++i)
  {
    convolution[i] = convolution_at_angular_velocity_point_0[i] + convolution_at_angular_velocity_point_1[i];
  } // j
  normalization = std::accumulate(projected_system_state.begin(), projected_system_state.end(), Real(0.0));
}