//
// Created by Nikita Kruk on 25.06.20.
//

#include "Convolution.hpp"
#include "FastDiscreteInverseFourierTransform.hpp"

#include <vector>
#include <iostream>
#include <algorithm>
#include <complex>

Convolution::Convolution() :
    f_transform_(kL, std::complex<Real>(0, 0)),
    g_transform_(kL, std::complex<Real>(0, 0)),
    convolution_transform_(kL, std::complex<Real>(0, 0))
{

}

Convolution::~Convolution()
{

}

void Convolution::SlowConvolution(const std::vector<Real> &f,
                                  const std::vector<Real> &g,
                                  std::vector<Real> &convolution)
{

}

void Convolution::FastConvolution(const std::vector<Real> &f,
                                  const std::vector<Real> &g,
                                  std::vector<Real> &convolution)
{
  FastDiscreteInverseFourierTransform fdift;
//  std::vector<std::complex<Real>> f_transform(kL, 0.0);
//  std::vector<std::complex<Real>> g_transform(kL, 0.0);
  fdift.MakeTransform(f, f_transform_);
  fdift.MakeTransform(g, g_transform_);

//  std::vector<std::complex<Real>> convolution_transform(kL, 0.0);
  for (int k = 0; k < kL; ++k)
  {
    convolution_transform_[k] = f_transform_[k] * g_transform_[k];
  } // k

  std::vector<std::complex<Real>> &convolution_double_transform = f_transform_;
//  std::fill(convolution_double_transform.begin(), convolution_double_transform.end(), std::complex<Real>(0.0, 0.0));
  fdift.MakeTransform(convolution_transform_, convolution_double_transform);
  for (int k = 0; k < kL; ++k)
  {
    convolution[utilities::PositiveModulo(-k, kL)] = convolution_double_transform[k].real() / kL;
  } // k
}