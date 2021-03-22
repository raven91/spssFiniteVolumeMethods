//
// Created by Nikita Kruk on 25.06.20.
//

#include "DynSysSplittingFastPtr.hpp"
#include "../FourierTransforms/Convolution.hpp"

#include <mpi.h>
#include <iostream>
#include <numeric> // std::accumulate
#include <algorithm> // std::fill

DynSysSplittingFastPtr::DynSysSplittingFastPtr(ThreadOmegaSharedMemory *thread) :
    DynSysSplittingPtr(thread),
    cos_kernel_at_spatial_point_(kL, 0.0),
    sin_kernel_at_spatial_point_(kL, 0.0)
{
  for (int i = 0; i < kL; ++i)
  {
    cos_kernel_at_spatial_point_[i] = (-kC2 + kC1) * std::cos(Phi(i) + kAlpha);
    sin_kernel_at_spatial_point_[i] = -kC1 * std::sin(Phi(i) + kAlpha);
  } // k
}

DynSysSplittingFastPtr::~DynSysSplittingFastPtr()
{

}

/*void DynSysSplittingFastPtr::CalculateFlux(Dimension dim)
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

  CalculateDensitySlopes(dim);
  VerifyPositivityOfDensityAtCellInterfaces(dim);

  MPI_Win &win_density_slope = thread_->GetWindow(std::string("density_slope_window"));
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win &win_flux = thread_->GetWindow(std::string("flux_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  MPI_Win_lock_all(0, win_flux);
  static std::vector<Real> convolution(kL, 0.0);
  Real normalization = 0.0;
  CalculateConvolution(convolution, normalization);
  for (const int &angular_velocity_index: thread_->GetAngularVelocityLoopIndices())
  {
    int j = angular_velocity_index;
    for (int i = 0; i < kL; ++i)
    {
      utilities::TwoDimIdxToOneDimIdx(i, j, idx);
      density_next_in_current_cell = rk_system_state_[idx] + delta / 2.0 * density_slope_[idx];
      switch (dim)
      {
        case Dimension::kPhi: utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
          break;
        case Dimension::kOmega: utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
          break;
      }
      density_prev_in_next_cell = rk_system_state_[idx_next] - delta / 2.0 * density_slope_[idx_next];

      Real velocity = 0.0;
      switch (dim)
      {
        case Dimension::kPhi: velocity = Omega(j);
          if (cfl_a_ < std::fabs(velocity))
          {
            cfl_a_ = std::fabs(velocity);
          }
          break;
        case Dimension::kOmega: CalculateVelocityAtCellInterface(i, j, velocity, convolution, normalization);
          if (cfl_b_ < std::fabs(velocity))
          {
            cfl_b_ = std::fabs(velocity);
          }
          break;
      }
      flux_[idx] = std::max((Real) 0.0, velocity) * density_next_in_current_cell
          + std::min((Real) 0.0, velocity) * density_prev_in_next_cell;
    } // k
  } // angular_velocity_index
//  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Win_sync(win_flux);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  MPI_Win_unlock_all(win_flux);
  thread_->SynchronizeVectorThroughoutClusters(flux_);
}*/

void DynSysSplittingFastPtr::CalculateAngularFlux()
{
  Real density_next_in_current_cell = 0.0, density_prev_in_next_cell = 0.0;
  int idx_next = 0, i = 0, j = 0;
  Real velocity = 0.0;

  CalculateDensitySlopes(Dimension::kPhi, std::string("density_slope_wrt_angle_window"), density_slope_wrt_angle_);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kPhi,
                                            std::string("density_slope_wrt_angle_window"),
                                            density_slope_wrt_angle_);
  thread_->SynchronizeVectorThroughoutClusters(density_slope_wrt_angle_);

  MPI_Win &win_density_slope = thread_->GetWindow(std::string("density_slope_wrt_angle_window"));
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win &win_flux = thread_->GetWindow(std::string("angular_flux_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  MPI_Win_lock_all(0, win_flux);
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    velocity = Omega(j);
    if (velocity > 0.0)
    {
      density_next_in_current_cell = rk_system_state_[idx] + kDphi / 2.0 * density_slope_wrt_angle_[idx];
      angular_flux_[idx] = velocity * density_next_in_current_cell;
    } else
    {
      utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
      density_prev_in_next_cell = rk_system_state_[idx_next] - kDphi / 2.0 * density_slope_wrt_angle_[idx_next];
      angular_flux_[idx] = velocity * density_prev_in_next_cell;
    }
    if (cfl_a_ < std::fabs(velocity))
    {
      cfl_a_ = std::fabs(velocity);
    }
  } // idx
//  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Win_sync(win_flux);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  MPI_Win_unlock_all(win_flux);
  thread_->SynchronizeVectorThroughoutClusters(angular_flux_);
}

void DynSysSplittingFastPtr::CalculateAngularVelocityFlux()
{
  Real density_next_in_current_cell = 0.0, density_prev_in_next_cell = 0.0;
  int idx_next = 0, i = 0, j = 0;
  Real velocity = 0.0;

  CalculateDensitySlopes(Dimension::kPhi, std::string("density_slope_wrt_angle_window"), density_slope_wrt_angle_);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kPhi,
                                            std::string("density_slope_wrt_angle_window"),
                                            density_slope_wrt_angle_);
  thread_->SynchronizeVectorThroughoutClusters(density_slope_wrt_angle_);
  CalculateDensitySlopes(Dimension::kOmega,
                         std::string("density_slope_wrt_angular_velocity_window"),
                         density_slope_wrt_angular_velocity_);
  VerifyPositivityOfDensityAtCellInterfaces(Dimension::kOmega,
                                            std::string("density_slope_wrt_angular_velocity_window"),
                                            density_slope_wrt_angular_velocity_);
  thread_->SynchronizeVectorThroughoutClusters(density_slope_wrt_angular_velocity_);

  MPI_Win &win_density_slope_wrt_angle = thread_->GetWindow(std::string("density_slope_wrt_angle_window"));
  MPI_Win &win_density_slope_wrt_angular_velocity =
      thread_->GetWindow(std::string("density_slope_wrt_angular_velocity_window"));
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win &win_flux = thread_->GetWindow(std::string("angular_velocity_flux_window"));
  MPI_Win_lock_all(0, win_density_slope_wrt_angle);
  MPI_Win_lock_all(0, win_density_slope_wrt_angular_velocity);
  MPI_Win_lock_all(0, win_rk_system_state);
  MPI_Win_lock_all(0, win_flux);
  static std::vector<Real> convolution(kL, 0.0);
  Real normalization = 0.0;
  CalculateConvolution(convolution, normalization);
  auto &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    CalculateVelocityAtCellInterface(i, j, velocity, convolution, normalization);
    if (velocity > 0.0)
    {
      density_next_in_current_cell = rk_system_state_[idx] + kDomega / 2.0 * density_slope_wrt_angular_velocity_[idx];
      angular_velocity_flux_[idx] = velocity * density_next_in_current_cell;
    } else
    {
      utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
      density_prev_in_next_cell =
          rk_system_state_[idx_next] - kDomega / 2.0 * density_slope_wrt_angular_velocity_[idx_next];
      angular_velocity_flux_[idx] = velocity * density_prev_in_next_cell;
    }
    if (cfl_b_ < std::fabs(velocity))
    {
      cfl_b_ = std::fabs(velocity);
    }
  } // idx
//  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Win_sync(win_flux);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope_wrt_angle);
  MPI_Win_unlock_all(win_density_slope_wrt_angular_velocity);
  MPI_Win_unlock_all(win_rk_system_state);
  MPI_Win_unlock_all(win_flux);
  thread_->SynchronizeVectorThroughoutClusters(angular_velocity_flux_);
}

/**
 * For the terms with convolution only
 * @param i
 * @param velocities
 */
void DynSysSplittingFastPtr::CalculateVelocityAtCellInterface(int i,
                                                              int j,
                                                              Real &velocity,
                                                              const std::vector<Real> &convolution,
                                                              Real normalization)
{
  Real potential_cur = 0.0, potential_next = 0.0;
  CalculateVelocityPotentialAtCellCenter(i, j, potential_cur, convolution, normalization);
  CalculateVelocityPotentialAtCellCenter(i,
                                         utilities::PositiveModulo(j + 1, kK),
                                         potential_next,
                                         convolution,
                                         normalization);
  velocity = -(potential_next - potential_cur) / kDomega;
}

/**
 * Calculation of a potential through convolution only
 * @param i
 * @param potentials
 */
void DynSysSplittingFastPtr::CalculateVelocityPotentialAtCellCenter(int i,
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
    Real log_density = std::log(rk_system_state_[idx]);
    if (std::isfinite(log_density))
    {
      potential += (kDiffusionConstant * log_density); // a fix for extremely small density values in diffusion terms
    }
  }
}

void DynSysSplittingFastPtr::CalculateConvolution(std::vector<Real> &convolution, Real &normalization)
{
  static std::vector<Real> projected_system_state(kL, 0.0), projected_density_slope(kL, 0.0);
  std::fill(projected_system_state.begin(), projected_system_state.end(), 0.0);
  std::fill(projected_density_slope.begin(), projected_density_slope.end(), 0.0);

  for (int m = 0; m < kK; ++m)
  {
    for (int l = 0; l < kL; ++l)
    {
      int idx = 0;
      utilities::TwoDimIdxToOneDimIdx(l, m, idx);
      projected_system_state[l] += rk_system_state_[idx];
      projected_density_slope[l] += density_slope_wrt_angle_[idx];
    } // l
  } // m

  // convolution of the interaction kernel
  Convolution convolution_engine;
  static std::vector<Real> convolution_at_angular_velocity_point_0(kL, 0.0);
  convolution_engine.FastConvolution(projected_system_state,
                                     sin_kernel_at_spatial_point_,
                                     convolution_at_angular_velocity_point_0);
  static std::vector<Real> convolution_at_angular_velocity_point_1(kL, 0.0);
  convolution_engine.FastConvolution(projected_density_slope,
                                     cos_kernel_at_spatial_point_,
                                     convolution_at_angular_velocity_point_1);

  for (int i = 0; i < kL; ++i)
  {
    convolution[i] = convolution_at_angular_velocity_point_0[i] + convolution_at_angular_velocity_point_1[i];
  } // j
  normalization = std::accumulate(projected_system_state.begin(), projected_system_state.end(), Real(0.0));
}