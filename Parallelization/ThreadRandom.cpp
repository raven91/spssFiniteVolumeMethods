//
// Created by Nikita Kruk on 29.06.20.
//

#include "ThreadRandom.hpp"
#include "Parallelization.hpp"

#include <numeric> // std::iota
#include <random>
#include <iostream>
#include <cassert>
#include <algorithm> // std::min_element, std::shuffle

ThreadRandom::ThreadRandom(int argc, char **argv) :
    ThreadOmegaSharedMemory(argc, argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
  assert(!(kK
      % number_of_mpich_threads_));  // width by height must be divisible by the number of threads for this kind of parallelization
  number_of_elements_per_mpich_thread_ = kL * kK / number_of_mpich_threads_;

  all_indices_ = std::vector<int>(kL * kK, 0);
  std::iota(all_indices_.begin(), all_indices_.end(), 0);
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(all_indices_.begin(), all_indices_.end(), g);
  loop_indices_ = std::vector<int>(number_of_elements_per_mpich_thread_, 0);
  std::copy(&all_indices_[rank_ * number_of_elements_per_mpich_thread_],
            &all_indices_[(rank_ + 1) * number_of_elements_per_mpich_thread_],
            loop_indices_.begin());
  // loop_indices_ is unrelated to angular_velocity_loop_indices_
  // angular_velocity_loop_indices_ is used for unique convolution computation
  std::vector<int> all_angular_velocity_indices(kK, 0);
  std::shuffle(all_angular_velocity_indices.begin(), all_angular_velocity_indices.end(), g);
  angular_velocity_loop_indices_ = std::vector<int>(kK / number_of_mpich_threads_, 0);
  std::copy(&all_angular_velocity_indices[rank_ * kK / number_of_mpich_threads_],
            &all_angular_velocity_indices[(rank_ + 1) * kK / number_of_mpich_threads_],
            angular_velocity_loop_indices_.begin());

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_communicator_);
  MPI_Comm_rank(shared_communicator_, &shared_rank_);
  MPI_Comm_size(shared_communicator_, &number_of_shared_mpich_threads_);
  shared_root_rank_ = 0;
  shared_ranks_ = std::vector<int>((unsigned long) number_of_mpich_threads_, 0);
  MPI_Gather(&shared_rank_, 1, MPI_INT, &shared_ranks_[0], 1, MPI_INT, root_rank_, MPI_COMM_WORLD);
  MPI_Bcast(&shared_ranks_[0], (int) shared_ranks_.size(), MPI_INT, root_rank_, MPI_COMM_WORLD);
  numbers_of_shared_mpich_threads_ = std::vector<int>((unsigned long) number_of_mpich_threads_, 0);
  MPI_Gather(&number_of_shared_mpich_threads_,
             1,
             MPI_INT,
             &numbers_of_shared_mpich_threads_[0],
             1,
             MPI_INT,
             root_rank_,
             MPI_COMM_WORLD);
  MPI_Bcast(&numbers_of_shared_mpich_threads_[0],
            (int) numbers_of_shared_mpich_threads_.size(),
            MPI_INT,
            root_rank_,
            MPI_COMM_WORLD);
  for (int i = 0; i < shared_ranks_.size(); ++i)
  {
    if (shared_ranks_[i] == shared_root_rank_)
    {
      shared_root_rank_indexes_.push_back(i);
    }
  } // i
}

ThreadRandom::~ThreadRandom()
{
  for (std::pair<const std::string, MPI_Win> &window : windows_)
  {
    MPI_Win_free(&(window.second));
  } // window
  MPI_Comm_free(&shared_communicator_);
}

int ThreadRandom::GetRootRank()
{
  return root_rank_;
}

int ThreadRandom::GetSharedRootRank()
{
  return shared_root_rank_;
}

const std::vector<int> &ThreadRandom::GetSharedRanks()
{
  return shared_ranks_;
}

const std::vector<int> &ThreadRandom::GetSharedRootRankIndexes()
{
  return shared_root_rank_indexes_;
}

const std::vector<int> &ThreadRandom::GetNumbersOfSharedMpichThreads()
{
  return numbers_of_shared_mpich_threads_;
}

MPI_Comm &ThreadRandom::GetSharedCommunicator()
{
  return shared_communicator_;
}

MPI_Win &ThreadRandom::GetWindow(const std::string &window_name)
{
  return windows_[window_name];
}

bool ThreadRandom::IsRoot()
{
  return (rank_ == root_rank_);
}

bool ThreadRandom::IsSharedRoot()
{
  return (shared_rank_ == shared_root_rank_);
}

void ThreadRandom::SynchronizeVector(std::vector<Real> &vec)
{
  Real *const p = &vec[0];
  SynchronizeVector(p, vec.size());
}

void ThreadRandom::SynchronizeVector(Real *const vec, long size)
{
  MPI_Win win;
  if (rank_ == root_rank_)
  {
    MPI_Win_create(&vec[0], size * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    for (int i : loop_indices_)
    {
      MPI_Put(&vec[i], 1, kRealTypeForMpi, root_rank_, i, 1, kRealTypeForMpi, win);
    } // i
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], size, kRealTypeForMpi, root_rank_, 0, size, kRealTypeForMpi, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadRandom::SynchronizeVectorThroughoutClusters(Real *const vec)
{
// if there are threads on different nodes
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(number_of_mpich_threads_ - 1);
      for (int i = 1; i < number_of_mpich_threads_; ++i)
      {
        std::vector<Real> buf(number_of_elements_per_mpich_thread_, 0.0);
        MPI_Irecv(&buf[0],
                  number_of_elements_per_mpich_thread_,
                  kRealTypeForMpi,
                  i,
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
        for (int j = 0; j < buf.size(); ++j)
        {
          int shuffled_index = all_indices_[number_of_elements_per_mpich_thread_ * i + j];
          vec[shuffled_index] = buf[j];
        } // j
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else
    {
      MPI_Request request;
      std::vector<Real> buf;
      for (int i : loop_indices_)
      {
        buf.push_back(vec[i]);
      } // i
      MPI_Isend(&buf[0],
                number_of_elements_per_mpich_thread_,
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    if (IsRoot())
    {
      std::vector<MPI_Request> requests(number_of_mpich_threads_ - 1);
      for (int i = 1; i < number_of_mpich_threads_; ++i)
      {
        std::vector<Real> buf(number_of_elements_per_mpich_thread_, 0.0);
        for (int j = 0; j < buf.size(); ++j)
        {
          int shuffled_index = all_indices_[number_of_elements_per_mpich_thread_ * i + j];
          buf[j] = vec[shuffled_index];
        } // j
        MPI_Isend(&buf[0],
                  number_of_elements_per_mpich_thread_,
                  kRealTypeForMpi,
                  i,
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else
    {
      MPI_Request request;
      std::vector<Real> buf(number_of_elements_per_mpich_thread_, 0.0);
      MPI_Irecv(&buf[0],
                number_of_elements_per_mpich_thread_,
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      for (int i = 0; i < number_of_elements_per_mpich_thread_; ++i)
      {
        int shuffled_index = loop_indices_[i];
        vec[shuffled_index] = buf[i];
      } // i
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadRandom::BroadcastVector(std::vector<Real> &vec)
{
  Real *const p = &vec[0];
  BroadcastVector(p, vec.size());
}

void ThreadRandom::BroadcastVector(Real *const vec, long size)
{
  MPI_Win win;
  if (IsRoot())
  {
    MPI_Win_create(&vec[0], size * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    for (int i : loop_indices_)
    {
      MPI_Get(&vec[i], 1, kRealTypeForMpi, root_rank_, 0, 1, kRealTypeForMpi, win);
    } // i
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadRandom::BroadcastValue(double &v)
{
  MPI_Bcast(&v, 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
}

void ThreadRandom::BroadcastVectorThroughoutClusters(Real *const vec)
{
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(number_of_mpich_threads_ - 1);
      for (int i = 1; i < number_of_mpich_threads_; ++i)
      {
        std::vector<Real> buf(number_of_elements_per_mpich_thread_, 0.0);
        for (int j = 0; j < buf.size(); ++j)
        {
          int shuffled_index = all_indices_[number_of_elements_per_mpich_thread_ * i + j];
          buf[j] = vec[shuffled_index];
        } // j
        MPI_Isend(&buf[0],
                  number_of_elements_per_mpich_thread_,
                  kRealTypeForMpi,
                  i,
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
      std::vector<Real> buf(number_of_elements_per_mpich_thread_, 0.0);
      MPI_Irecv(&buf[0],
                number_of_elements_per_mpich_thread_,
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      for (int i = 0; i < number_of_elements_per_mpich_thread_; ++i)
      {
        int shuffled_index = loop_indices_[i];
        vec[shuffled_index] = buf[i];
      } // i
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadRandom::FindMinValue(double &v)
{
  std::vector<Real> vec((unsigned long) (number_of_mpich_threads_), 0.0);
  MPI_Gather(&v, 1, kRealTypeForMpi, &vec[0], 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
  if (IsRoot())
  {
    v = *std::min_element(vec.begin(), vec.end());
  }
  MPI_Bcast(&v, 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
}

void ThreadRandom::ComputeSumIntoRootOnly(std::vector<int> &vec)
{
  std::vector<int> global_vec(vec.size(), 0);
  MPI_Reduce(&vec[0], &global_vec[0], (int) vec.size(), MPI_INT, MPI_SUM, root_rank_, MPI_COMM_WORLD);
  if (IsRoot())
  {
    vec.swap(global_vec);
  }
}

void ThreadRandom::BroadcastCondition(int &condition)
{
  MPI_Bcast(&condition, 1, MPI_INT, root_rank_, MPI_COMM_WORLD);
}

void ThreadRandom::AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name)
{
  if (windows_.find(window_name) == windows_.end())
  {
    MPI_Win win;
    MPI_Win_allocate_shared(size * sizeof(Real),
                            sizeof(Real),
                            MPI_INFO_NULL,
                            shared_communicator_,
                            &vec,
                            &win);
    windows_[window_name] = win;
  } else
  {
    std::cout << "Window '" << window_name << "' already exists. The previous window is kept." << std::endl;
  }
}

void ThreadRandom::FreeSharedWindow(const std::string &window_name)
{
  if (windows_.find(window_name) != windows_.end())
  {
    MPI_Win_free(&(windows_[window_name]));
    windows_.erase(window_name);
  } else
  {
    std::cout << "Window '" << window_name << "' does not exist. Nothing is erased." << std::endl;
  }
}