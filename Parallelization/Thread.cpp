//
// Created by Nikita Kruk on 24.06.20.
//

#include "Thread.hpp"

#include <cassert>
#include <iostream>
#include <numeric>  // std::iota
#include <algorithm> // std::min_element

Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  rank_ = 0;
  number_of_mpich_threads_ = 1;
  assert(!(kK
      % number_of_mpich_threads_));  // omega grid must be divisible by the number of threads for this kind of parallelization
  number_of_elements_per_mpich_thread_ = kL * kK / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(kL * kK / number_of_mpich_threads_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_elements_per_mpich_thread_);
  angular_velocity_loop_indices_ = std::vector<int>(kK / number_of_mpich_threads_, 0);
  std::iota(angular_velocity_loop_indices_.begin(),
            angular_velocity_loop_indices_.end(),
            rank_ * kK / number_of_mpich_threads_);
}

Thread::~Thread()
{

}

int Thread::GetRank()
{
  return rank_;
}

int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}

int Thread::GetNumberOfElementsPerMpichThread()
{
  return number_of_elements_per_mpich_thread_;
}

const std::vector<int> &Thread::GetLoopIndices()
{
  return loop_indices_;
}

const std::vector<int> &Thread::GetAngularVelocityLoopIndices()
{
  return angular_velocity_loop_indices_;
}

bool Thread::IsRoot()
{
  return true;
}

void Thread::SynchronizeVector(std::vector<Real> &vec)
{

}

void Thread::BroadcastVector(std::vector<Real> &vec)
{

}

void Thread::BroadcastValue(double &v)
{

}

void Thread::FindMinValue(double &v)
{

}

void Thread::ComputeSumIntoRootOnly(std::vector<int> &vec)
{

}

void Thread::BroadcastCondition(int &condition)
{

}

//void Thread::InitializeNeighborThreads(std::vector<std::vector<int>> &spatial_neighbor_indices)
//{
//
//}