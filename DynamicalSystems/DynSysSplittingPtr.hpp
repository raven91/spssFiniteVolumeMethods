//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadOmegaSharedMemory.hpp"

#include <vector>

/**
 * Dynamical system computes the RHS of a given PDE
 * using the Strang splitting for
 * angular update (\Delta t / 2), angular velocity update (\Delta t), and angular update (\Delta t / 2) steps
 */
class DynSysSplittingPtr
{
 public:

  explicit DynSysSplittingPtr(ThreadOmegaSharedMemory *thread);
  virtual ~DynSysSplittingPtr();

  Real CalculateAngularUpdate(const Real *const system_state,
                              const std::vector<Real> &k_prev,
                              std::vector<Real> &k_next,
                              Real k_coef,
                              Real dt);
  Real CalculateAngularVelocityUpdate(const Real *const system_state,
                                      const std::vector<Real> &k_prev,
                                      std::vector<Real> &k_next,
                                      Real k_coef,
                                      Real dt);

  void ResetFluxLimiterCount();
  const std::vector<int> &GetFluxLimiterCount() const
  { return flux_limiter_count_; }
  Real GetCflA() const
  { return cfl_a_; }
  Real GetCflB() const
  { return cfl_b_; }

  ThreadOmegaSharedMemory *thread_;
  std::vector<Real> sin_of_phase_;
  std::vector<Real> cos_of_phase_;
  Real *rk_system_state_;
//  Real *density_slope_;
  Real *density_slope_wrt_angle_;
  Real *density_slope_wrt_angular_velocity_;
//  Real *flux_;
  Real *angular_flux_;
  Real *angular_velocity_flux_;
  long vector_size_;
  Real cfl_a_;
  Real cfl_b_;
  std::vector<int> flux_limiter_count_;

//  virtual void CalculateFlux(Dimension dim) = 0;
  virtual void CalculateAngularFlux() = 0;
  virtual void CalculateAngularVelocityFlux() = 0;
  void CalculateDensitySlopes(Dimension dim,
                              const std::string &density_slope_window_name,
                              Real *const density_slope);
  void CalculateDensitySlopesWrtAngle();
  void CalculateDensitySlopesWrtAngularVelocity();
  void VerifyPositivityOfDensityAtCellInterfaces(Dimension dim,
                                                 const std::string &density_slope_window_name,
                                                 Real *const density_slope);
  void RecalculateDensitySlopeWithFluxLimiter(Dimension dim,
                                              int idx);

 private:

};

#endif //SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGPTR_HPP
