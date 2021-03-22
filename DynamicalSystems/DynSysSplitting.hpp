//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTING_HPP
#define SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTING_HPP

#include "../Definitions.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

/**
 * Dynamical system computes the RHS of a given PDE
 * using the Strang splitting for
 * angular update (\Delta t / 2), angular velocity update (\Delta t), and angular update (\Delta t / 2) steps
 */
class DynSysSplitting
{
 public:

  explicit DynSysSplitting(Thread *thread);
  virtual ~DynSysSplitting();

  Real CalculateAngularUpdate(const std::vector<Real> &system_state,
                              const std::vector<Real> &k_prev,
                              std::vector<Real> &k_next,
                              Real k_coef,
                              Real dt);
  Real CalculateAngularVelocityUpdate(const std::vector<Real> &system_state,
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

  Thread *thread_;
  std::vector<Real> sin_of_phase_;
  std::vector<Real> cos_of_phase_;
//  std::vector<Real> density_slope_;
  std::vector<Real> density_slope_wrt_angle_;
  std::vector<Real> density_slope_wrt_angular_velocity_;
  Real cfl_a_;
  Real cfl_b_;
  std::vector<int> flux_limiter_count_;

//  virtual void CalculateFlux(Dimension dim, const std::vector<Real> &system_state, std::vector<Real> &flux) = 0;
  virtual void CalculateAngularFlux(const std::vector<Real> &system_state, std::vector<Real> &angular_flux) = 0;
  virtual void CalculateAngularVelocityFlux(const std::vector<Real> &system_state,
                                            std::vector<Real> &angular_velocity_flux) = 0;
  void CalculateDensitySlopes(Dimension dim, const std::vector<Real> &system_state);
  void VerifyPositivityOfDensityAtCellInterfaces(Dimension dim,
                                                 const std::vector<Real> &system_state,
                                                 const std::vector<Real> &density_slope);
  void RecalculateDensitySlopeWithFluxLimiter(Dimension dim,
                                              const std::vector<Real> &system_state,
                                              int idx);

 private:

};

#endif //SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTING_HPP
