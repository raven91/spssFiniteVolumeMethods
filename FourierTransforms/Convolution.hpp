//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_CONVOLUTION_HPP
#define SPSSFINITEVOLUMEMETHODS_CONVOLUTION_HPP

#include "../Definitions.hpp"

#include <vector>
#include <complex>

class Convolution
{
 public:

  Convolution();
  ~Convolution();

  void SlowConvolution(const std::vector<Real> &f, const std::vector<Real> &g, std::vector<Real> &convolution);
  void FastConvolution(const std::vector<Real> &f, const std::vector<Real> &g, std::vector<Real> &convolution);

 private:

  std::vector<std::complex<Real>> f_transform_;
  std::vector<std::complex<Real>> g_transform_;
  std::vector<std::complex<Real>> convolution_transform_;

};

#endif //SPSSFINITEVOLUMEMETHODS_CONVOLUTION_HPP
