//
// Created by Nikita Kruk on 25.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_BINARYOBSERVERPTR_HPP
#define SPSSFINITEVOLUMEMETHODS_BINARYOBSERVERPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadOmegaSharedMemory.hpp"

#include <string>
#include <fstream>
#include <chrono>

class BinaryObserverPtr
{
 public:

  int should_terminate_;

  explicit BinaryObserverPtr(ThreadOmegaSharedMemory *thread);
  ~BinaryObserverPtr();

  void SaveSystemState(Real *const system_state, long size, Real t);
  void SaveSummaryStatistics(const Real *const system_state, long size, Real t);
  void SaveAdditionalInformation(const std::vector<int> &flux_limiter_count,
                                 Real cfl_a,
                                 Real cfl_b,
                                 Real t);

 protected:

  ThreadOmegaSharedMemory *thread_;

 private:

  std::string output_file_name_;
  std::ofstream output_file_;
  std::string summary_statistics_file_name_;
  std::ofstream summary_statistics_file_;
  std::string flux_limiter_activity_file_name_;
  std::ofstream flux_limiter_activity_file_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int output_time_counter_[3];
  int output_time_threshold_[3];

  void ValidateSolution(Real *const system_state, long size, Real t);

};

#endif //SPSSFINITEVOLUMEMETHODS_BINARYOBSERVERPTR_HPP
