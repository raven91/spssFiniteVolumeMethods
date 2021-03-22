//
// Created by Nikita Kruk on 24.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_THREADOMEGASHAREDMEMORY_HPP
#define SPSSFINITEVOLUMEMETHODS_THREADOMEGASHAREDMEMORY_HPP

#include "Thread.hpp"

#include <mpi.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

class ThreadOmegaSharedMemory : public Thread
{
 public:

  ThreadOmegaSharedMemory(int argc, char **argv);
  ~ThreadOmegaSharedMemory();

  virtual int GetRootRank();
  virtual int GetSharedRootRank();
  virtual const std::vector<int> &GetSharedRanks();
  virtual const std::vector<int> &GetSharedRootRankIndexes();
  virtual const std::vector<int> &GetNumbersOfSharedMpichThreads();
  virtual MPI_Comm &GetSharedCommunicator();
  virtual MPI_Win &GetWindow(const std::string &window_name);

  virtual bool IsRoot();
  virtual bool IsSharedRoot();

  void SynchronizeVector(Dimension dim, std::vector<Real> &vec);
  virtual void SynchronizeVector(std::vector<Real> &vec);
  virtual void SynchronizeVector(Real *const vec, long size);
  virtual void SynchronizeVectorThroughoutClusters(Real *const vec);
  virtual void BroadcastVector(std::vector<Real> &vec);
  virtual void BroadcastVector(Real *const vec, long size);
  virtual void BroadcastVectorThroughoutClusters(Real *const vec);
  virtual void BroadcastValue(double &v);
  virtual void FindMinValue(double &v);
  virtual void ComputeSumIntoRootOnly(std::vector<int> &vec);
  virtual void BroadcastCondition(int &condition);
//  void InitializeNeighborThreads(std::vector<std::vector<int>> &spatial_neighbor_indices) override;

  virtual void AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name);
//  void AllocateSharedWindow(MPI_Aint size, Real *vec, MPI_Win &window);
  virtual void FreeSharedWindow(const std::string &window_name);
//  void FreeSharedWindow(MPI_Win &window);

 protected:

  int shared_rank_;
  std::vector<int> shared_ranks_;
  int shared_root_rank_;
  std::vector<int> shared_root_rank_indexes_;
  int number_of_shared_mpich_threads_; // in a shared memory communicator
  std::vector<int> numbers_of_shared_mpich_threads_; // in each shared memory communicator
  MPI_Comm shared_communicator_;
  std::unordered_map<std::string, MPI_Win> windows_;
//  std::unordered_set<int> neighbor_threads_;
  std::unordered_map<int, int> index_to_prev_omega_thread_;

};

#endif //SPSSFINITEVOLUMEMETHODS_THREADOMEGASHAREDMEMORY_HPP
