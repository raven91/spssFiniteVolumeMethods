//
// Created by Nikita Kruk on 25.06.20.
//

#include "DynSysSplittingPtr.hpp"

#include <mpi.h>
#include <algorithm> // std::min
#include <cassert>

DynSysSplittingPtr::DynSysSplittingPtr(ThreadOmegaSharedMemory *thread) :
    thread_(thread),
    sin_of_phase_(kL, 0.0),
    cos_of_phase_(kL, 0.0),
    vector_size_(kL * kK),
    flux_limiter_count_(kDim, 0)
{
  for (int k = 0; k < kL; ++k)
  {
    sin_of_phase_[k] = std::sin(Phi(k));
    cos_of_phase_[k] = std::cos(Phi(k));
  } // k

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(vector_size_, rk_system_state_, std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread->AllocateSharedWindow(0, rk_system_state_, std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("rk_system_state_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &rk_system_state_);
    assert(vector_size_ == (size / sizeof(Real)));
  }
//  MPI_Barrier(thread_->GetSharedCommunicator());

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(vector_size_,
                                  density_slope_wrt_angle_,
                                  std::string("density_slope_wrt_angle_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread->AllocateSharedWindow(0, density_slope_wrt_angle_, std::string("density_slope_wrt_angle_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("density_slope_wrt_angle_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &density_slope_wrt_angle_);
    assert(vector_size_ == (size / sizeof(Real)));
  }
//  MPI_Barrier(thread_->GetSharedCommunicator());

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(vector_size_,
                                  density_slope_wrt_angular_velocity_,
                                  std::string("density_slope_wrt_angular_velocity_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread->AllocateSharedWindow(0,
                                 density_slope_wrt_angular_velocity_,
                                 std::string("density_slope_wrt_angular_velocity_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("density_slope_wrt_angular_velocity_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &density_slope_wrt_angular_velocity_);
    assert(vector_size_ == (size / sizeof(Real)));
  }
//  MPI_Barrier(thread_->GetSharedCommunicator());

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(vector_size_, angular_flux_, std::string("angular_flux_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread->AllocateSharedWindow(0, angular_flux_, std::string("angular_flux_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("angular_flux_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &angular_flux_);
    assert(vector_size_ == (size / sizeof(Real)));
  }
//  MPI_Barrier(thread_->GetSharedCommunicator());

  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(vector_size_, angular_velocity_flux_, std::string("angular_velocity_flux_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread->AllocateSharedWindow(0, angular_velocity_flux_, std::string("angular_velocity_flux_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("angular_velocity_flux_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &angular_velocity_flux_);
    assert(vector_size_ == (size / sizeof(Real)));
  }
//  MPI_Barrier(thread_->GetSharedCommunicator());
}

DynSysSplittingPtr::~DynSysSplittingPtr()
{
  sin_of_phase_.clear();
  cos_of_phase_.clear();

  MPI_Barrier(thread_->GetSharedCommunicator());
  thread_->FreeSharedWindow(std::string("rk_system_state_window"));
  thread_->FreeSharedWindow(std::string("density_slope_wrt_angle_window"));
  thread_->FreeSharedWindow(std::string("density_slope_wrt_angular_velocity_window"));
  thread_->FreeSharedWindow(std::string("angular_flux_window"));
  thread_->FreeSharedWindow(std::string("angular_velocity_flux_window"));
}

Real DynSysSplittingPtr::CalculateAngularUpdate(const Real *const system_state,
                                                const std::vector<Real> &k_prev,
                                                std::vector<Real> &k_next,
                                                Real k_coef,
                                                Real dt)
{
  cfl_a_ = 0.0;

  /*MPI_Win &win_sys = thread_->GetWindow(std::string("system_state_window"));
  MPI_Win_lock_all(0, win_sys);*/
  MPI_Win &win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[idx] = system_state[idx] + k_coef * dt * k_prev[idx];
  } // idx
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  /*MPI_Win_unlock_all(win_sys);*/
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

//  CalculateFlux(Dimension::kPhi);
  CalculateAngularFlux();

  int i = 0, j = 0;
  int idx_prev_phi = 0;
  MPI_Win &win_flux = thread_->GetWindow(std::string("angular_flux_window"));
  MPI_Win_lock_all(0, win_flux);
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev_phi);
    k_next[idx] = -(angular_flux_[idx] - angular_flux_[idx_prev_phi]) / kDphi;
  } // idx
  MPI_Win_unlock_all(win_flux);

  /*if (thread_->IsRoot())
  {
    std::cout << "CFL: " << dt << " <= " << std::min(kDx / (4.0 * cfl_a_), kDy / (4.0 * cfl_b_)) << std::endl;
  }*/
  return kDphi / (2.0 * cfl_a_);
}

Real DynSysSplittingPtr::CalculateAngularVelocityUpdate(const Real *const system_state,
                                                        const std::vector<Real> &k_prev,
                                                        std::vector<Real> &k_next,
                                                        Real k_coef,
                                                        Real dt)
{
  cfl_b_ = 0.0;

  /*MPI_Win &win_sys = thread_->GetWindow(std::string("system_state_window"));
  MPI_Win_lock_all(0, win_sys);*/
  MPI_Win &win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    rk_system_state_[idx] = system_state[idx] + k_coef * dt * k_prev[idx];
  } // idx
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  /*MPI_Win_unlock_all(win_sys);*/
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

//  std::fill(flux_.begin(), flux_.end(), 0.0);
//  CalculateFlux(Dimension::kOmega);
  CalculateAngularVelocityFlux();

  int i = 0, j = 0;
  int idx_prev_omega = 0;
  MPI_Win &win_flux = thread_->GetWindow(std::string("angular_velocity_flux_window"));
  MPI_Win_lock_all(0, win_flux);
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev_omega);
    // impose zero flux boundary conditions
    if (j == 0)
    {
      k_next[idx] = -(angular_velocity_flux_[idx] - 0.0) / kDomega;
    } else if (j == kK - 1)
    {
      k_next[idx] = -(0.0 - angular_velocity_flux_[idx_prev_omega]) / kDomega;
    } else
    {
      k_next[idx] = -(angular_velocity_flux_[idx] - angular_velocity_flux_[idx_prev_omega]) / kDomega;
    }
  } // idx
  MPI_Win_unlock_all(win_flux);

  /*if (thread_->IsRoot())
  {
    std::cout << "CFL: " << dt << " <= " << kDphi / (2.0 * cfl_c_) << std::endl;
  }*/
  return kDphi / (2.0 * cfl_b_);
}

//void DynSysSplittingPtr::CalculateFlux(Dimension dim)
//{
//
//}

void DynSysSplittingPtr::CalculateAngularFlux()
{

}

void DynSysSplittingPtr::CalculateAngularVelocityFlux()
{

}

/**
 *
 * @param dim : 0, 1
 * @param system_state
 */
void DynSysSplittingPtr::CalculateDensitySlopes(Dimension dim,
                                                const std::string &density_slope_window_name,
                                                Real *const density_slope)
{
  int i = 0, j = 0;
  int idx_next = 0, idx_prev = 0;
  Real delta = 0.0;
  MPI_Win &win_density_slope = thread_->GetWindow(density_slope_window_name);
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    switch (dim)
    {
      case Dimension::kPhi : utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
        utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev);
        delta = kDphi;
        break;
      case Dimension::kOmega : utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
        utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev);
        delta = kDomega;
        break;
    }
    density_slope[idx] = (rk_system_state_[idx_next] - rk_system_state_[idx_prev]) / (2.0 * delta);
  } // idx
  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  thread_->SynchronizeVectorThroughoutClusters(density_slope);
}

void DynSysSplittingPtr::CalculateDensitySlopesWrtAngle()
{
  int i = 0, j = 0;
  int idx_next = 0, idx_prev = 0;
  MPI_Win &win_density_slope = thread_->GetWindow(std::string("density_slope_wrt_angle_window"));
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
    utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev);
    density_slope_wrt_angle_[idx] = (rk_system_state_[idx_next] - rk_system_state_[idx_prev]) / (2.0 * kDphi);
  } // idx
  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  thread_->SynchronizeVectorThroughoutClusters(density_slope_wrt_angle_);
}

void DynSysSplittingPtr::CalculateDensitySlopesWrtAngularVelocity()
{
  int i = 0, j = 0;
  int idx_next = 0, idx_prev = 0;
  MPI_Win &win_density_slope = thread_->GetWindow(std::string("density_slope_wrt_angular_velocity_window"));
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);
    utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
    utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev);
    density_slope_wrt_angular_velocity_[idx] =
        (rk_system_state_[idx_next] - rk_system_state_[idx_prev]) / (2.0 * kDomega);
  } // idx
  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  thread_->SynchronizeVectorThroughoutClusters(density_slope_wrt_angular_velocity_);
}

/**
 *
 * @param dim : 0, 1
 * @param system_state
 */
void DynSysSplittingPtr::VerifyPositivityOfDensityAtCellInterfaces(Dimension dim,
                                                                   const std::string &density_slope_window_name,
                                                                   Real *const density_slope)
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
  MPI_Win &win_density_slope = thread_->GetWindow(density_slope_window_name);
  MPI_Win &win_rk_system_state = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_density_slope);
  MPI_Win_lock_all(0, win_rk_system_state);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int idx : loop_indices)
  {
    density_next = rk_system_state_[idx] + delta / 2.0 * density_slope[idx];
    density_prev = rk_system_state_[idx] - delta / 2.0 * density_slope[idx];
    if ((density_prev < 0.0) || (density_next < 0.0))
    {
      RecalculateDensitySlopeWithFluxLimiter(dim, idx);
    }
//    else if (density_next < 0.0)
//    {
//      RecalculateDensitySlopeWithFluxLimiter(dim, idx, false);
//    }
  } // idx
  MPI_Win_sync(win_density_slope);
//  MPI_Win_sync(win_rk_system_state);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_density_slope);
  MPI_Win_unlock_all(win_rk_system_state);
  thread_->SynchronizeVectorThroughoutClusters(density_slope);
}

void DynSysSplittingPtr::RecalculateDensitySlopeWithFluxLimiter(Dimension dim, int idx)
{
  int i = 0, j = 0;
  int idx_next = 0, idx_prev = 0;
  utilities::OneDimIdxToTwoDimIdx(idx, i, j);

  Real theta = 1.0;
  switch (dim)
  {
    case Dimension::kPhi : utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i + 1, kL), j, idx_next);
      utilities::TwoDimIdxToOneDimIdx(utilities::PositiveModulo(i - 1, kL), j, idx_prev);
      density_slope_wrt_angle_[idx] =
          utilities::Minmod(theta * (rk_system_state_[idx_next] - rk_system_state_[idx]) / kDphi,
                            density_slope_wrt_angle_[idx],
                            theta * (rk_system_state_[idx] - rk_system_state_[idx_prev]) / kDphi);
      break;
    case Dimension::kOmega : utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j + 1, kK), idx_next);
      utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev);
      density_slope_wrt_angular_velocity_[idx] =
          utilities::Minmod(theta * (rk_system_state_[idx_next] - rk_system_state_[idx]) / kDomega,
                            density_slope_wrt_angular_velocity_[idx],
                            theta * (rk_system_state_[idx] - rk_system_state_[idx_prev]) / kDomega);
      break;
  }

//  density_slope_[idx] = utilities::Minmod(theta * (rk_system_state_[idx_next] - rk_system_state_[idx]) / delta,
//                                          density_slope_[idx],
//                                          theta * (rk_system_state_[idx] - rk_system_state_[idx_prev]) / delta);
  /*// max possible density slope
  if (left_is_negative)
  {
    density_slope_[idx] = rk_system_state_[idx] / (delta / 2.0);
  } else
  {
    density_slope_[idx] = -rk_system_state_[idx] / (delta / 2.0);
  }*/

  switch (dim)
  {
    case Dimension::kPhi : ++flux_limiter_count_[0];
      break;
    case Dimension::kOmega : ++flux_limiter_count_[1];
      break;
  }
}

void DynSysSplittingPtr::ResetFluxLimiterCount()
{
  std::fill(flux_limiter_count_.begin(), flux_limiter_count_.end(), 0);
}