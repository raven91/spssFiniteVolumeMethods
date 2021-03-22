//
// Created by Nikita Kruk on 29.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_THREADRANDOM_HPP
#define SPSSFINITEVOLUMEMETHODS_THREADRANDOM_HPP

#include "ThreadOmegaSharedMemory.hpp"

#include <mpi.h>
#include <string>
#include <unordered_map>

class ThreadRandom : public ThreadOmegaSharedMemory
{
 public:

  ThreadRandom(int argc, char **argv);
  ~ThreadRandom();

  virtual int GetRootRank();
  virtual int GetSharedRootRank();
  virtual const std::vector<int> &GetSharedRanks();
  virtual const std::vector<int> &GetSharedRootRankIndexes();
  virtual const std::vector<int> &GetNumbersOfSharedMpichThreads();
  virtual MPI_Comm &GetSharedCommunicator();
  virtual MPI_Win &GetWindow(const std::string &window_name);

  virtual bool IsRoot();
  virtual bool IsSharedRoot();

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

  void AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name);
  void FreeSharedWindow(const std::string &window_name);

 private:

  std::vector<int> all_indices_;

};

#endif //SPSSFINITEVOLUMEMETHODS_THREADRANDOM_HPP
