//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_INITIALCONDITIONS_HPP
#define SPSSFINITEVOLUMEMETHODS_INITIALCONDITIONS_HPP

#include "../Definitions.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <cassert>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>

#if defined(MULTIPRECISION)
#include <boost/multiprecision/mpfr.hpp>
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>> MultiprecisionReal;
#else
typedef Real MultiprecisionReal;
#endif

class InitialConditions
{
 public:

  template<typename T>
  static void UniformPerturbed(T &vec, long size)
  {
    std::mt19937 mersenne_twister_generator(std::random_device{}());
    std::uniform_real_distribution<Real> unif_real_dist(0, 0.0001);
//    boost::math::normal_distribution<Real> norm_dist(0.0, std::sqrt(kDiffusionConstant / kFriction));

    int i = 0, j = 0;
    for (int idx = 0; idx < size; ++idx)
    {
      utilities::OneDimIdxToTwoDimIdx(idx, i, j);
      vec[idx] = 1.0 / (kPhiSize) * (1.0 / kOmegaSize) + unif_real_dist(mersenne_twister_generator);
//      vec[idx] =
//          (1 / kPhiSize) * boost::math::pdf(norm_dist, Omega(j)) + unif_real_dist(mersenne_twister_generator);
    } // idx
    std::cout << typeid(vec).name() << std::endl;
  }

  template<typename T>
  static void VonMisesGaussianProduct(T &vec, long size)
  {
#if defined(MULTIPRECISION)
    MultiprecisionReal::default_precision(50);
    using namespace boost::multiprecision;
#else
    using namespace std;
#endif

    const Real op_magnitude = 0.99496192622318713;
    const MultiprecisionReal gamma = kSigma * kFriction * op_magnitude / kDiffusionConstant;
    const Real two_pi = boost::math::constants::two_pi<Real>();
    const MultiprecisionReal s = std::sqrt(kDiffusionConstant / kFriction);

    int idx = 0;
    for (int i = 0; i < kL; ++i)
    {
      for (int j = 0; j < kK; ++j)
      {
        utilities::TwoDimIdxToOneDimIdx(i, j, idx);
        Real phi_i = Phi(i);
        Real omega_j = Omega(j);
        MultiprecisionReal density_value = exp(gamma * std::cos(phi_i)) / (two_pi * boost::math::cyl_bessel_i(0, gamma))
            / sqrt(two_pi * s * s) * exp(-0.5 * omega_j * omega_j / (s * s));
        vec[idx] = (Real) density_value;
      } // j
    } // i
  }

  template<typename T>
  static void UniformGaussianProductPerturbed(T &vec, long size)
  {
    std::mt19937 mersenne_twister_generator(std::random_device{}());
    std::uniform_real_distribution<Real> unif_real_dist(0, 0.0001);
    const Real two_pi = boost::math::constants::two_pi<Real>();
    const Real s = std::sqrt(kDiffusionConstant / kFriction);

    int idx = 0;
    for (int i = 0; i < kL; ++i)
    {
      for (int j = 0; j < kK; ++j)
      {
        utilities::TwoDimIdxToOneDimIdx(i, j, idx);
        Real omega_j = Omega(j);
        vec[idx] = 1.0 / (std::pow(two_pi, 1.5) * s) * std::exp(-0.5 * omega_j * omega_j / (s * s))
            * (1.0 + unif_real_dist(mersenne_twister_generator));
      } // j
    } // i
  }

  template<typename T>
  static void CoarseGrainedParticleDensity(T &vec, long size)
  {
    std::string folder("/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/");
    std::ifstream particle_dynamics_file
        (folder + "v0_1_xi_0.1_sigma_1_rho_0.8_alpha_0.3_Dphi_0_N_100000_0_0.bin", std::ios::binary | std::ios::in);
    assert(particle_dynamics_file.is_open());
    const int n = 100000, s = 4;
    const int skip_time = 0;
    particle_dynamics_file.seekg(skip_time * (1l + s * n) * sizeof(float), std::ios::beg);
    std::vector<float> particle_solution(s * n, 0.0f);
    float time = 0.0f;
    particle_dynamics_file.read((char *) &time, sizeof(float));
    particle_dynamics_file.read((char *) &particle_solution[0], s * n * sizeof(float));
    particle_dynamics_file.close();

    for (int i = 0; i < n; ++i)
    {
      Real phi_i = particle_solution[s * i + 2];
      phi_i -= std::floor(phi_i * kPhiRSize) * kPhiSize;
      Real omega_i = particle_solution[s * i + 3];
//      omega_i -= std::floor(omega_i * kOmegaRSize) * kOmegaSize;

      int cell_phi = int(phi_i / kPhiSize * kL);
      int cell_omega = int(omega_i / kOmegaSize * kK + kK / 2.0);
      int cell_idx = 0;
      utilities::TwoDimIdxToOneDimIdx(cell_phi, cell_omega, cell_idx);
      vec[cell_idx] += 1.0;
    } // i
  }

};

#endif //SPSSFINITEVOLUMEMETHODS_INITIALCONDITIONS_HPP
