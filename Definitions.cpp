//
// Created by Nikita Kruk on 24.06.20.
//

#include "Definitions.hpp"

#include <vector>

namespace utilities
{
  int PositiveModulo(int i, int n)
  {
    if (i < 0)
    {
      return n + i;
    } else if (i >= n)
    {
      return i - n;
    } else
    {
      return i;
    }
//    return (i % n + n) % n;
//  return (i + n) % n;
  }

  void TwoDimIdxToOneDimIdx(int phi, int omega, int &idx)
  {
    // the winding order is phi->omega
    idx = phi + kL * omega;
  }

  void TwoDimIdxToOneDimIdx(int phi, int omega, int &idx, int num_cells_phi, int num_cells_omega)
  {
    // the winding order is phi->omega
    idx = phi + num_cells_phi * omega;
  }

  void OneDimIdxToTwoDimIdx(int idx, int &phi, int &omega)
  {
    // the winding order is phi->omega
    omega = idx / kL;
    phi = idx % kL;
  }

  void OneDimIdxToTwoDimIdx(int idx, int &phi, int &omega, int num_cells_phi, int num_cells_omega)
  {
    // the winding order is phi->x->y
    omega = idx / num_cells_phi;
    phi = idx % num_cells_phi;
  }

  Real Minmod(const Real &a, const Real &b, const Real &c)
  {
    if (a > 0.0 && b > 0.0 && c > 0.0)
    {
      return std::min(std::min(a, b), c);
    } else if (a < 0.0 && b < 0.0 && c < 0.0)
    {
      return std::max(std::max(a, b), c);
    } else
    {
      return Real(0.0);
    }
  }
} // namespace utilities

/**
 *
 * @param i
 * @return x_i is the grid value in the middle of a cell, x_{i-0.5} and x_{i+0.5} are the values on the cell interface
 * 0 = x_{-0.5} \equiv x_{N-0.5} = 1
 */
Real Phi(const Real &i)
{
  return i * kDphi;
}

Real Omega(const Real &j)
{
  return (-kK / 2 + j) * kDomega;
}