//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFASTPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFASTPTR_HPP

#include "DynSysSplittingPtr.hpp"
#include "../Parallelization/ThreadOmegaSharedMemory.hpp"

#include <vector>

class DynSysSplittingFastPtr : public DynSysSplittingPtr
{
 public:

  explicit DynSysSplittingFastPtr(ThreadOmegaSharedMemory *thread);
  ~DynSysSplittingFastPtr();

 private:

  std::vector<Real> cos_kernel_at_spatial_point_;
  std::vector<Real> sin_kernel_at_spatial_point_;

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

#endif //SPSSFINITEVOLUMEMETHODS_DYNSYSSPLITTINGFASTPTR_HPP
