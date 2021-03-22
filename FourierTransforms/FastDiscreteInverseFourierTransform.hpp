//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_FASTDISCRETEINVERSEFOURIERTRANSFORM_HPP
#define SPSSFINITEVOLUMEMETHODS_FASTDISCRETEINVERSEFOURIERTRANSFORM_HPP

#include "../Definitions.hpp"

#include <vector>
#include <complex>
#include <bitset>

/**
 * Implementation based on:
 * Knuth, D. The Art of Computer Programming
 */
class FastDiscreteInverseFourierTransform
{
 public:

  FastDiscreteInverseFourierTransform();
  ~FastDiscreteInverseFourierTransform();

  static void MakeTransform(const std::vector<Real> &original, std::vector<std::complex<Real>> &transform);
  static void MakeTransform(const std::vector<std::complex<Real>> &original, std::vector<std::complex<Real>> &transform);

 private:

  // \omega^0, \omega^1, ... , \omega^{2^k-1}
  static std::vector<std::complex<Real>> omega_;
  static std::vector<int> reversed_indexes_;

  template<std::size_t N>
  void Reverse(std::bitset<N> &b);

};

#endif //SPSSFINITEVOLUMEMETHODS_FASTDISCRETEINVERSEFOURIERTRANSFORM_HPP
