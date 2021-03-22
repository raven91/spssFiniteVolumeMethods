//
// Created by Nikita Kruk on 25.06.20.
//

#include "BinaryObserverPtr.hpp"

#include <mpi.h>
#include <sstream>
#include <cassert>
#include <iostream>
#include <numeric>  // std::accumulate
#include <iterator> // std::copy
#include <complex>

BinaryObserverPtr::BinaryObserverPtr(ThreadOmegaSharedMemory *thread) :
    thread_(thread),
    output_time_counter_{0, 0, 0},
    output_time_threshold_{1 * kInverseDt, kInverseDt, kInverseDt}, // mod 1 - save at every dt
    should_terminate_(0)
{
  if (thread->IsRoot())
  {
    integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
    std::string folder("/work/scratch/nk59zoce/cpp/spssFiniteVolumeMethods/");
#elif defined(__linux__) && defined(BCS_CLUSTER)
    std::string folder("/home/nkruk/cpp/spssFiniteVolumeMethods/output/");
#elif defined(__linux__)
    std::string folder("/home/nikita/Documents/spssFiniteVolumeMethods/");
#elif defined(__APPLE__)
//    std::string folder("/Users/nikita/Documents/Projects/spss/spssFiniteVolumeMethods/");
    std::string folder("/Volumes/Kruk/spss/spssFiniteVolumeMethods/continuation/prolonged_dynamics/");
#endif
    std::ostringstream output_file_name_buffer;
    output_file_name_buffer << folder << "dt_" << kDt << "_v0_" << kMicroscopicVelocity << "_xi_" << kFriction
                            << "_sigma_" << kSigma << "_rho_" << kRho << "_alpha_" << kAlpha
                            << "_Dphi_" << kDiffusionConstant << "_" << kL << "_" << kK << ".bin";
    output_file_name_ = output_file_name_buffer.str();
    if (thread->IsRoot())
    {
      std::remove(output_file_name_.c_str());
    }
    output_file_.open(output_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(output_file_.is_open());

    std::ostringstream summary_statistics_file_name_buffer;
    summary_statistics_file_name_buffer << folder << "dt_" << kDt << "_v0_" << kMicroscopicVelocity << "_xi_"
                                        << kFriction
                                        << "_sigma_" << kSigma << "_rho_" << kRho << "_alpha_" << kAlpha
                                        << "_Dphi_" << kDiffusionConstant << "_" << kL << "_" << kK << ".txt";
    summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
    if (thread->IsRoot())
    {
      std::remove(summary_statistics_file_name_.c_str());
    }
    summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
    assert(summary_statistics_file_.is_open());

    std::ostringstream flux_limiter_activity_file_name_buffer;
    flux_limiter_activity_file_name_buffer << folder << "flux_limiter_dt_" << kDt << "_v0_" << kMicroscopicVelocity
                                           << "_xi_" << kFriction
                                           << "_sigma_" << kSigma << "_rho_" << kRho << "_alpha_" << kAlpha
                                           << "_Dphi_" << kDiffusionConstant << "_" << kL << "_" << kK << ".txt";
    flux_limiter_activity_file_name_ = flux_limiter_activity_file_name_buffer.str();
    if (thread->IsRoot())
    {
      std::remove(flux_limiter_activity_file_name_.c_str());
    }
    flux_limiter_activity_file_.open(flux_limiter_activity_file_name_, std::ios::out | std::ios::app);
    assert(flux_limiter_activity_file_.is_open());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

BinaryObserverPtr::~BinaryObserverPtr()
{
  if (thread_->IsRoot())
  {
    if (output_file_.is_open())
    {
      output_file_.close();
    }
    if (summary_statistics_file_.is_open())
    {
      summary_statistics_file_.close();
    }
    if (flux_limiter_activity_file_.is_open())
    {
      flux_limiter_activity_file_.close();
    }
  }
}

void BinaryObserverPtr::SaveSystemState(Real *const system_state, long size, Real t)
{
  thread_->SynchronizeVector(system_state, size);
  if (thread_->IsRoot())
  {
    ValidateSolution(system_state, size, t);
    std::cout << "total mass: "
              << std::accumulate(system_state, system_state + size, (Real) 0.0) * kDphi * kDomega << std::endl;

    std::chrono::duration<RealForOutput> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "t:" << t << " | integration time:" << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

    if (!(output_time_counter_[0] % output_time_threshold_[0]))
    {
      RealForOutput t_float = RealForOutput(t);
      static std::vector<RealForOutput> system_state_float(size, 0.0);
      std::copy(&system_state[0], &system_state[size], &system_state_float[0]);
      output_file_.write((char *) &t_float, sizeof(RealForOutput));
      output_file_.write((char *) &system_state_float[0], system_state_float.size() * sizeof(RealForOutput));
    }
    ++output_time_counter_[0];
  }
}

void BinaryObserverPtr::ValidateSolution(Real *const system_state, long size, Real t)
{
  for (int idx = 0; idx < size; ++idx)
  {
    if (!std::isfinite(system_state[idx]))
    {
      std::cout << "finiteness assertion failed!" << std::endl;
      should_terminate_ = 1;
      return;
    }

    if (system_state[idx] < 0.0)
    {
      std::cout << "positivity assertion failed: density == " << system_state[idx] << " < 0.0!" << std::endl;
      should_terminate_ = 1;
      return;
    }
  } // idx
}

void BinaryObserverPtr::SaveSummaryStatistics(const Real *const system_state, long size, Real t)
{
  if (thread_->IsRoot())
  {
    if (!(output_time_counter_[1] % output_time_threshold_[1]))
    {
      float t_float = (float) t;
      std::complex<float> polar_order_parameter(0.0f, 0.0f);
      std::complex<float> nematic_order_parameter(0.0f, 0.0f);
      float normalization = 0.0;
      int i = 0, j = 0;
      for (int idx = 0; idx < size; ++idx)
      {
        utilities::OneDimIdxToTwoDimIdx(idx, i, j);
        polar_order_parameter +=
            std::complex<float>(std::cos((float) Phi(i)), std::sin((float) Phi(i))) * float(system_state[idx]);
        nematic_order_parameter += std::complex<float>(std::cos(2 * float(Phi(i))),
                                                       std::sin(2 * float(Phi(i)))) * float(system_state[idx]);
        normalization += system_state[idx];
      } // i
      polar_order_parameter /= normalization;
      nematic_order_parameter /= normalization;

      summary_statistics_file_ << t_float << '\t' << std::abs(polar_order_parameter) << '\t'
                               << std::arg(polar_order_parameter) << '\t'
                               << std::abs(nematic_order_parameter) << '\t'
                               << std::arg(nematic_order_parameter);
      summary_statistics_file_ << std::endl;
    }
    ++output_time_counter_[1];
  }
}

// number of applications of flux limiters
// velocities used for CFL conditions
void BinaryObserverPtr::SaveAdditionalInformation(const std::vector<int> &flux_limiter_count,
                                                  Real cfl_a,
                                                  Real cfl_b,
                                                  Real t)
{
  if (thread_->IsRoot())
  {
    if (!(output_time_counter_[2] % output_time_threshold_[2]))
    {
      flux_limiter_activity_file_ << t << '\t';
      for (int c : flux_limiter_count)
      {
        flux_limiter_activity_file_ << c << '\t';
      } // fl
      flux_limiter_activity_file_ << cfl_a << '\t' << cfl_b;
      flux_limiter_activity_file_ << std::endl;
    }
    ++output_time_counter_[2];
  }
}