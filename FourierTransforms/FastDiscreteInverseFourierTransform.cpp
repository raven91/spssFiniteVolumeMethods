//
// Created by Nikita Kruk on 25.06.20.
//

#include "FastDiscreteInverseFourierTransform.hpp"

#include <bitset>
#include <iostream>

std::vector<std::complex<Real>> FastDiscreteInverseFourierTransform::omega_(0);
std::vector<int> FastDiscreteInverseFourierTransform::reversed_indexes_(0);

FastDiscreteInverseFourierTransform::FastDiscreteInverseFourierTransform()
{
  if (omega_.empty())
  {
    std::vector<std::complex<Real>> basic_omega(kBits);
    basic_omega[0] = -1.0; // \omega_1
    basic_omega[1] = std::complex<Real>(0.0, 1.0); // \omega_2
    for (int r = 2; r < kBits; ++r)
    {
      basic_omega[r].real(std::sqrt((1.0 + basic_omega[r - 1].real()) / 2.0));
      basic_omega[r].imag(basic_omega[r - 1].imag() / (2.0 * basic_omega[r].real()));
    } // r

    omega_ = std::vector<std::complex<Real>>(kL);
    for (int j = 0; j < kL; ++j)
    {
      omega_[j] = 1.0;
      std::bitset<kBits> j_binary(j);
      for (int i = 0; i < j_binary.size(); ++i)
      {
        if (j_binary[i])
        {
          omega_[j] *= basic_omega[kBits - i - 1];
        }
      } // i
    } // j
  }

  if (reversed_indexes_.empty())
  {
    reversed_indexes_ = std::vector<int>(kL);
    for (int idx = 0; idx < kL; ++idx)
    {
      std::bitset<kBits> idx_binary(idx);
      Reverse(idx_binary);
      reversed_indexes_[idx] = (int) idx_binary.to_ulong();
    } // idx
  }
}

FastDiscreteInverseFourierTransform::~FastDiscreteInverseFourierTransform()
{

}

void FastDiscreteInverseFourierTransform::MakeTransform(const std::vector<Real> &original,
                                                        std::vector<std::complex<Real>> &transform)
{
  static std::vector<std::complex<Real>> original_complex(original.size());
  for (int i = 0; i < original.size(); ++i)
  {
    original_complex[i] = original[i];
  } // i
  MakeTransform(original_complex, transform);
}

void FastDiscreteInverseFourierTransform::MakeTransform(const std::vector<std::complex<Real>> &original,
                                                        std::vector<std::complex<Real>> &transform)
{
  static std::vector<std::complex<Real>> A_prev(kL), A_curr(kL);

  // Pass 0
  for (int idx = 0; idx < kL; ++idx)
  {
    A_curr[idx] = original[idx];
  } // t

  // Pass 1 to k
  int idx_0 = 0, idx_1 = 0, idx_omega = 0;
  for (int pass = 1; pass <= kBits; ++pass)
  {
    A_prev = A_curr;
    for (int idx = 0; idx < kL; ++idx)
    {
      std::bitset<kBits> idx_binary(idx);
      idx_binary[kBits - pass] = 0;
      idx_0 = (int) idx_binary.to_ulong();
      idx_binary[kBits - pass] = 1;
      idx_1 = (int) idx_binary.to_ulong();

//      std::bitset<kBits> idx_binary_omega(idx);
//      idx_binary_omega >>= (kBits - pass);
//      Reverse(idx_binary_omega);
//      idx_omega = (int) idx_binary_omega.to_ulong();
      idx_omega = reversed_indexes_[idx >> (kBits - pass)];
      A_curr[idx] = A_prev[idx_0] + omega_[idx_omega] * A_prev[idx_1];
    } // idx
  } // pass

  for (int idx = 0; idx < kL; ++idx)
  {
//    std::bitset<kBits> idx_binary(idx);
//    Reverse(idx_binary);
//    transform[idx_binary.to_ulong()] = A_curr[idx];
    transform[reversed_indexes_[idx]] = A_curr[idx];
  } // idx
}

template<std::size_t N>
void FastDiscreteInverseFourierTransform::Reverse(std::bitset<N> &b)
{
  bool t;
  for (std::size_t i = 0; i < N / 2; ++i)
  {
    t = b[i];
    b[i] = b[N - i - 1];
    b[N - i - 1] = t;
  }
}