//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFAST_HPP
#define SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFAST_HPP

#include "DynSysSplitting.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class DynSysSplittingFast : public DynSysSplitting
{
 public:

  explicit DynSysSplittingFast(Thread *thread);
  ~DynSysSplittingFast();

 private:

  std::vector<Real> cos_kernel_at_spatial_point_;
  std::vector<Real> sin_kernel_at_spatial_point_;

//  virtual void CalculateFlux(Dimension dim, const std::vector<Real> &system_state, std::vector<Real> &flux);
  virtual void CalculateAngularFlux(const std::vector<Real> &system_state, std::vector<Real> &angular_flux);
  virtual void CalculateAngularVelocityFlux(const std::vector<Real> &system_state,
                                            std::vector<Real> &angular_velocity_flux);
  void CalculateVelocityAtCellInterface(const std::vector<Real> &system_state,
                                        int i,
                                        int j,
                                        Real &velocity,
                                        const std::vector<Real> &convolution,
                                        Real normalization);
  void CalculateVelocityPotentialAtCellCenter(const std::vector<Real> &system_state,
                                              int i,
                                              int j,
                                              Real &potential,
                                              const std::vector<Real> &convolution,
                                              Real normalization);
  void CalculateConvolution(const std::vector<Real> &system_state, std::vector<Real> &convolution, Real &normalization);

};

#endif //SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFAST_HPP
