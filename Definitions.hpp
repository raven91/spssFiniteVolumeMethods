//
// Created by Nikita Kruk on 24.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_DEFINITIONS_HPP
#define SPSSFINITEVOLUMEMETHODS_DEFINITIONS_HPP

#define MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY
//#define MPI_PARAMETER_SCAN
//#define MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM

//#define LICHTENBERG
//#define BCS_CLUSTER

#include <cmath>
#include <vector>

typedef double Real; // should be consistent with MPI_DATATYPE
typedef double RealForOutput;

const int kDim = 2; // angle + angular velocity
enum class Dimension
{
  kPhi, kOmega
};

// model-specific constants
const Real kPhiSize = 2.0 * M_PI;
const Real kOmegaSize = 20.0;
const Real kPhiRSize = 1.0 / kPhiSize;
const Real kOmegaRSize = 1.0 / kOmegaSize;

const Real kMicroscopicVelocity = 1.0;
const Real kFriction = 0.1;
const Real kSigma = 1.0;
const Real kRho = 0.8;
const Real kAlpha = 0.3;
//extern Real kAlpha;
const Real kDiffusionConstant = 0.01;
//extern Real kDiffusionConstant;

// discretization-specific constants
const int kBits = 7;
const int kL = (int) std::pow(2.0, kBits);
const int kK = 400;

const Real kDphi = kPhiSize / kL;
const Real kDomega = kOmegaSize / kK;
const Real kT0 = 0.0;
const Real kT1 = 1000.0;
const Real kDt = 0.005;
const int kInverseDt = 1.0 / kDt;

// implementation-specific constants
const Real kC0 = kDphi / 2.0;
const Real kC1 = std::sin(kDphi / 2.0) / (kDphi / 2.0);
const Real kC2 = std::cos(kDphi / 2.0);
const Real kC3 = std::sin(kDphi / 2.0);

namespace utilities
{
  int PositiveModulo(int i, int n);
  void TwoDimIdxToOneDimIdx(int phi, int omega, int &idx);
  void TwoDimIdxToOneDimIdx(int phi, int omega, int &idx, int num_cells_phi, int num_cells_omega);
  void OneDimIdxToTwoDimIdx(int idx, int &phi, int &omega);
  void OneDimIdxToTwoDimIdx(int idx, int &phi, int &omega, int num_cells_phi, int num_cells_omega);
  Real Minmod(const Real &a, const Real &b, const Real &c);
} // namespace utilities
Real Phi(const Real &i);
Real Omega(const Real &j);

#endif //SPSSFINITEVOLUMEMETHODS_DEFINITIONS_HPP
