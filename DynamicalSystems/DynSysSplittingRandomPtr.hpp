//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGRANDOMPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGRANDOMPTR_HPP

#include "DynSysSplittingPtr.hpp"
#include "../Parallelization/ThreadRandom.hpp"

#include <vector>

class DynSysSplittingRandomPtr : public DynSysSplittingPtr
{
 public:

  explicit DynSysSplittingRandomPtr(ThreadRandom *thread);
  ~DynSysSplittingRandomPtr();

 private:

  std::vector<Real> cos_kernel_at_spatial_point_;
  std::vector<Real> sin_kernel_at_spatial_point_;
  std::vector<std::vector<Real>> convolutions_;
  std::vector<Real> normalizations_;

//  virtual void CalculateFlux(Dimension dim);
  virtual void CalculateAngularFlux();
  virtual void CalculateAngularVelocityFlux();
  void CalculateVelocityAtCellInterface(int i,
                                        int j,
                                        Real &velocity,
                                        const std::vector<Real> &convolution,
                                        Real normalization);
  void CalculateVelocityPotentialAtCellCenter(int i,
                                              int j,
                                              Real &potential,
                                              const std::vector<Real> &convolution,
                                              Real normalization);
  void CalculateConvolution(std::vector<Real> &convolution, Real &normalization);
};

#endif //SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGRANDOMPTR_HPP
