//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_BINARYOBSERVER_HPP
#define SPSSFINITEVOLUMEMETHODS_BINARYOBSERVER_HPP

#include "../Definitions.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>
#include <fstream>
#include <chrono>

class BinaryObserver
{
 public:

  int should_terminate_;

  explicit BinaryObserver(Thread *thread);
  ~BinaryObserver();

  void SaveSystemState(const std::vector<Real> &system_state, Real t);
  void SaveSummaryStatistics(const std::vector<Real> &system_state, Real t);
  void SaveAdditionalInformation(const std::vector<int> &flux_limiter_count,
                                 Real cfl_a,
                                 Real cfl_b,
                                 Real t);

 protected:

  Thread *thread_;

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

  void ValidateSolution(const std::vector<Real> &system_state);

};

#endif //SPSSFINITEVOLUMEMETHODS_BINARYOBSERVER_HPP
