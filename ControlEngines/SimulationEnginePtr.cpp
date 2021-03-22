//
// Created by Nikita Kruk on 25.06.20.
//

#include "SimulationEnginePtr.hpp"
#include "InitialConditions.hpp"
#include "../Observers/BinaryObserverPtr.hpp"
#include "../Steppers/RungeKutta2StepperWithSplittingPtr.hpp"
#include "../DynamicalSystems/DynSysSplittingFastPtr.hpp"

#include <algorithm> // std::copy
#include <random>

SimulationEnginePtr::SimulationEnginePtr(ThreadOmegaSharedMemory *thread) :
    thread_(thread),
    system_state_size_(kL * kK)
{
  system_state_ = new Real[system_state_size_];
}

SimulationEnginePtr::~SimulationEnginePtr()
{
  delete[] system_state_;
}

void SimulationEnginePtr::InitializeSystemStateByRule()
{
  if (thread_->IsRoot())
  {
    InitialConditions::UniformPerturbed<Real *const>(system_state_, system_state_size_);
//    InitialConditions::VonMisesGaussianProduct<Real *const>(system_state_, system_state_size_);
//    InitialConditions::UniformGaussianProductPerturbed<Real *const>(system_state_, system_state_size_);
//    InitialConditions::CoarseGrainedParticleDensity<Real *const>(system_state_, system_state_size_);

    // normalize the probability density function
    Real overall_density =
        std::accumulate(&system_state_[0], &system_state_[system_state_size_], Real(0.0)) * kDphi * kDomega;
    std::for_each(&system_state_[0], &system_state_[system_state_size_], [&](Real &ss) { ss /= overall_density; });
  }

  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEnginePtr::InitializeSystemStateFromFile(bool should_add_perturbation)
{
  if (thread_->IsRoot())
  {
#if defined(__linux__) && defined(BCS_CLUSTER)
    std::ifstream init_cond_file("/home/nkruk/cpp/spssFiniteVolumeMethods/input/dt_0.005_v0_1_xi_0.1_sigma_1_rho_0.8_alpha_0.3_Dphi_0.01_128_400.bin",
         std::ios::in | std::ios::binary);
#elif defined(__APPLE__)
    std::ifstream init_cond_file
        ("/Volumes/Kruk/spss/spssFiniteVolumeMethods/continuation/prolonged_dynamics/dt_0.005_v0_1_xi_0.1_sigma_1_rho_0.8_alpha_0.3_Dphi_0.0125_128_400.bin",
         std::ios::in | std::ios::binary);
#endif
    assert(init_cond_file.is_open());
    int t_init = 1000, dt_recip = 1;
    init_cond_file.seekg(t_init * dt_recip * (1l + system_state_size_) * sizeof(RealForOutput), std::ios::cur);
    RealForOutput t = 0.0;
    init_cond_file.read((char *) &t, sizeof(RealForOutput));
    std::vector<RealForOutput> input_system_state(system_state_size_, RealForOutput(0.0));
    init_cond_file.read((char *) &input_system_state[0], system_state_size_ * sizeof(RealForOutput));
    init_cond_file.close();
    std::copy(input_system_state.begin(), input_system_state.end(), &system_state_[0]);
    if (should_add_perturbation)
    {
      AddPerturbationToSystemState();
    }
  }
  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEnginePtr::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

//  InitializeSystemStateByRule();
  InitializeSystemStateFromFile(true);
  BinaryObserverPtr observer(thread_);
  DynSysSplittingFastPtr system(thread_);
  RungeKutta2StepperWithSplittingPtr stepper(thread_, kDt);

  Real t = kT0;
  observer.SaveSystemState(system_state_, system_state_size_, t);
  while (t <= kT1)
  {
    t += stepper.GetDt();
    stepper.DoStep(system, system_state_);
    observer.SaveSystemState(system_state_, system_state_size_, t);
    observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    observer.SaveAdditionalInformation(stepper.GetAverageFluxLimiterCount(),
                                       system.GetCflA(),
                                       system.GetCflB(),
                                       t);
    thread_->BroadcastCondition(observer.should_terminate_);
    if (observer.should_terminate_)
    {
      break;
    }
  }

  if (thread_->IsRoot())
  {
    std::cout << "simulation ended" << std::endl;
  }
}

void SimulationEnginePtr::RunContinuationMethod()
{
  /*if (thread_->IsRoot())
  {
    std::cout << "continuation method started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads"
              << std::endl;
  }

  Real initial_diffusion = 0.01;
  InitializeSystemStateFromFile();
  for (int i_diffusion = 0; i_diffusion < 17; ++i_diffusion)
  {
    kDiffusionConstant = initial_diffusion + 0.005 * i_diffusion;
    AddPerturbationToSystemState();
    BinaryObserverPtr observer(thread_);
    DynSysSplittingFastPtr system(thread_);
    RungeKutta2StepperWithSplittingPtr stepper(thread_, kDt);

    Real t = kT0;
    observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= kT1)
    {
      t += stepper.GetDt();
      stepper.DoStep(system, system_state_);
      observer.SaveSystemState(system_state_, system_state_size_, t);
      observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
      observer.SaveAdditionalInformation(stepper.GetAverageFluxLimiterCount(),
                                         system.GetCflA(),
                                         system.GetCflB(),
                                         t);
      thread_->BroadcastCondition(observer.should_terminate_);
      if (observer.should_terminate_)
      {
        break;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    thread_->BroadcastVector(system_state_, system_state_size_);
  } // i_diffusion

  if (thread_->IsRoot())
  {
    std::cout << "continuation method ended" << std::endl;
  }*/
}

void SimulationEnginePtr::AddPerturbationToSystemState()
{
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::uniform_real_distribution<Real> unif_real_dist(0, 0.0001);
  for (int idx = 0; idx < system_state_size_; ++idx)
  {
    system_state_[idx] += unif_real_dist(mersenne_twister_generator);
  } // idx
  // normalize the probability density function
  Real overall_density =
      std::accumulate(&system_state_[0], &system_state_[system_state_size_], Real(0.0)) * kDphi * kDomega;
  std::for_each(&system_state_[0], &system_state_[system_state_size_], [&](Real &ss) { ss /= overall_density; });
}