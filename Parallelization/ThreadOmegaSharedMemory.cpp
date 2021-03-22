//
// Created by Nikita Kruk on 24.06.20.
//

#include "ThreadOmegaSharedMemory.hpp"
#include "Parallelization.hpp"

#include <mpi.h>
#include <numeric>
#include <cassert>
#include <algorithm> // std::min_element
#include <utility>
#include <iostream>

ThreadOmegaSharedMemory::ThreadOmegaSharedMemory(int argc, char **argv) :
    Thread(argc, argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
  assert(!(kK
      % number_of_mpich_threads_));  // width by height must be divisible by the number of threads for this kind of parallelization
  number_of_elements_per_mpich_thread_ = kL * kK / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(number_of_elements_per_mpich_thread_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_elements_per_mpich_thread_);
  angular_velocity_loop_indices_ = std::vector<int>(kK / number_of_mpich_threads_, 0);
  std::iota(angular_velocity_loop_indices_.begin(),
            angular_velocity_loop_indices_.end(),
            rank_ * kK / number_of_mpich_threads_);

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

ThreadOmegaSharedMemory::~ThreadOmegaSharedMemory()
{
  for (std::pair<const std::string, MPI_Win> &window : windows_)
  {
    MPI_Win_free(&(window.second));
  } // window
  MPI_Comm_free(&shared_communicator_);
}

int ThreadOmegaSharedMemory::GetRootRank()
{
  return root_rank_;
}

const std::vector<int> &ThreadOmegaSharedMemory::GetSharedRanks()
{
  return shared_ranks_;
}

const std::vector<int> &ThreadOmegaSharedMemory::GetSharedRootRankIndexes()
{
  return shared_root_rank_indexes_;
}

const std::vector<int> &ThreadOmegaSharedMemory::GetNumbersOfSharedMpichThreads()
{
  return numbers_of_shared_mpich_threads_;
}

int ThreadOmegaSharedMemory::GetSharedRootRank()
{
  return shared_root_rank_;
}

MPI_Comm &ThreadOmegaSharedMemory::GetSharedCommunicator()
{
  return shared_communicator_;
}

MPI_Win &ThreadOmegaSharedMemory::GetWindow(const std::string &window_name)
{
  return windows_[window_name];
}

bool ThreadOmegaSharedMemory::IsRoot()
{
  return (rank_ == root_rank_);
}

bool ThreadOmegaSharedMemory::IsSharedRoot()
{
  return (shared_rank_ == shared_root_rank_);
}

void ThreadOmegaSharedMemory::SynchronizeVector(Dimension dim, std::vector<Real> &vec)
{
  /*MPI_Win win;
  std::vector<Real> vec_prev
      (&vec[rank_ * number_of_elements_per_mpich_thread_], &vec[(rank_ + 1) * number_of_elements_per_mpich_thread_]);
  MPI_Win_create(&vec_prev[0], vec_prev.size() * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  MPI_Win_fence(0, win);
  if (Dimension::kOmega == dim)
  {
    for (const std::pair<const int, int> &global_idx_thrd : index_to_prev_omega_thread_)
    {
      int global_idx = 0 + kL * global_idx_thrd.first;
      int local_idx = global_idx - global_idx_thrd.second * number_of_elements_per_mpich_thread_;
      MPI_Get(&vec[global_idx],
              kL,
              kRealTypeForMpi,
              global_idx_thrd.second,
              local_idx,
              kL,
              kRealTypeForMpi,
              win);
    } // idx_thrd
  }
  MPI_Win_fence(0, win);
  MPI_Win_free(&win);*/
}

void ThreadOmegaSharedMemory::SynchronizeVector(std::vector<Real> &vec)
{
  Real *const p = &vec[0];
  SynchronizeVector(p, vec.size());
}

void ThreadOmegaSharedMemory::SynchronizeVector(Real *const vec, long size)
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
    MPI_Put(&vec[rank_ * number_of_elements_per_mpich_thread_],
            number_of_elements_per_mpich_thread_,
            kRealTypeForMpi,
            root_rank_,
            rank_ * number_of_elements_per_mpich_thread_,
            number_of_elements_per_mpich_thread_,
            kRealTypeForMpi,
            win);
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], size, kRealTypeForMpi, root_rank_, 0, size, kRealTypeForMpi, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadOmegaSharedMemory::SynchronizeVectorThroughoutClusters(Real *const vec)
{
  // if there are threads on different nodes
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Recv(&vec[number_of_elements_per_mpich_thread_ * shared_root_rank_indexes_[i]],
//                 number_of_elements_per_mpich_thread_ * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
//                 kRealTypeForMpi,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD,
//                 MPI_STATUS_IGNORE);
        MPI_Irecv(&vec[number_of_elements_per_mpich_thread_ * shared_root_rank_indexes_[i]],
                  number_of_elements_per_mpich_thread_ * numbers_of_shared_mpich_threads_[shared_root_rank_indexes_[i]],
                  kRealTypeForMpi,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Send(&vec[number_of_elements_per_mpich_thread_ * rank_],
//               number_of_elements_per_mpich_thread_ * numbers_of_shared_mpich_threads_[rank_],
//               kRealTypeForMpi,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD);
      MPI_Isend(&vec[number_of_elements_per_mpich_thread_ * rank_],
                number_of_elements_per_mpich_thread_ * numbers_of_shared_mpich_threads_[rank_],
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
//        MPI_Send(&vec[0],
//                 number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
//                 kRealTypeForMpi,
//                 shared_root_rank_indexes_[i],
//                 0,
//                 MPI_COMM_WORLD);
        MPI_Isend(&vec[0],
                  number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
                  kRealTypeForMpi,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
//      MPI_Recv(&vec[0],
//               number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
//               kRealTypeForMpi,
//               root_rank_,
//               0,
//               MPI_COMM_WORLD,
//               MPI_STATUS_IGNORE);
      MPI_Irecv(&vec[0],
                number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadOmegaSharedMemory::BroadcastVector(std::vector<Real> &vec)
{
  MPI_Win win;
  if (rank_ == root_rank_)
  {
    MPI_Win_create(&vec[0], vec.size() * sizeof(Real), sizeof(Real), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Win_fence(0, win);
  } else
  {
    MPI_Win_create(nullptr, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Get(&vec[0], vec.size(), kRealTypeForMpi, root_rank_, 0, vec.size(), kRealTypeForMpi, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadOmegaSharedMemory::BroadcastVector(Real *const vec, long size)
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
    MPI_Get(&vec[0], size, kRealTypeForMpi, root_rank_, 0, size, kRealTypeForMpi, win);
    MPI_Win_fence(0, win);
  }
  MPI_Win_free(&win);
}

void ThreadOmegaSharedMemory::BroadcastValue(double &v)
{
  MPI_Bcast(&v, 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
}

void ThreadOmegaSharedMemory::BroadcastVectorThroughoutClusters(Real *const vec)
{
  if (shared_root_rank_indexes_.size() > 1)
  {
    if (IsRoot())
    {
      std::vector<MPI_Request> requests(shared_root_rank_indexes_.size() - 1);
      for (int i = 1; i < shared_root_rank_indexes_.size(); ++i)
      {
        MPI_Isend(&vec[0],
                  number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
                  kRealTypeForMpi,
                  shared_root_rank_indexes_[i],
                  0,
                  MPI_COMM_WORLD,
                  &requests[i - 1]);
      } // i
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE);
    } else if (IsSharedRoot())
    {
      MPI_Request request;
      MPI_Irecv(&vec[0],
                number_of_elements_per_mpich_thread_ * number_of_mpich_threads_,
                kRealTypeForMpi,
                root_rank_,
                0,
                MPI_COMM_WORLD,
                &request);
      MPI_Waitall(1, &request, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD); // for all the non-root processes
  }
}

void ThreadOmegaSharedMemory::FindMinValue(double &v)
{
  std::vector<Real> vec((unsigned long) (number_of_mpich_threads_), 0.0);
  MPI_Gather(&v, 1, kRealTypeForMpi, &vec[0], 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
  if (IsRoot())
  {
    v = *std::min_element(vec.begin(), vec.end());
  }
  MPI_Bcast(&v, 1, kRealTypeForMpi, root_rank_, MPI_COMM_WORLD);
}

void ThreadOmegaSharedMemory::ComputeSumIntoRootOnly(std::vector<int> &vec)
{
  std::vector<int> global_vec(vec.size(), 0);
  MPI_Reduce(&vec[0], &global_vec[0], (int) vec.size(), MPI_INT, MPI_SUM, root_rank_, MPI_COMM_WORLD);
  if (IsRoot())
  {
    vec.swap(global_vec);
  }
}

void ThreadOmegaSharedMemory::BroadcastCondition(int &condition)
{
  MPI_Bcast(&condition, 1, MPI_INT, root_rank_, MPI_COMM_WORLD);
}

/*// parameter is not used
void ThreadOmegaSharedMemory::InitializeNeighborThreads(std::vector<std::vector<int>> &spatial_neighbor_indices)
{
  int i = 0, j = 0;
  int idx_prev_omega = 0;
  int rank_prev_omega = 0;
  for (const int &idx : loop_indices_)
  {
    utilities::OneDimIdxToTwoDimIdx(idx, i, j);

    utilities::TwoDimIdxToOneDimIdx(i, utilities::PositiveModulo(j - 1, kK), idx_prev_omega);
    rank_prev_omega = idx_prev_omega / number_of_elements_per_mpich_thread_;
    if (rank_prev_omega != rank_)
    {
      index_to_prev_omega_thread_[i + kL * utilities::PositiveModulo(j - 1, kK)] = rank_prev_omega;
    }
  } // idx
}*/

void ThreadOmegaSharedMemory::AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name)
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

void ThreadOmegaSharedMemory::FreeSharedWindow(const std::string &window_name)
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