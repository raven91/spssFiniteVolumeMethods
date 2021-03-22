//
// Created by Nikita Kruk on 24.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_THREAD_HPP
#define SPSSFINITEVOLUMEMETHODS_THREAD_HPP

#include "../Definitions.hpp"

#include <vector>

class Thread
{
 public:

  Thread(int argc, char **argv);
  virtual ~Thread();

  int GetRank();
  int GetNumberOfMpichThreads();
  int GetNumberOfElementsPerMpichThread();
  const std::vector<int> &GetLoopIndices();
  const std::vector<int> &GetAngularVelocityLoopIndices();

  virtual bool IsRoot();

  virtual void SynchronizeVector(std::vector<Real> &vec);
  virtual void BroadcastVector(std::vector<Real> &vec);
  virtual void BroadcastValue(double &v);
  virtual void FindMinValue(double &v);
  virtual void ComputeSumIntoRootOnly(std::vector<int> &vec);
  virtual void BroadcastCondition(int &condition);

//  virtual void InitializeNeighborThreads(std::vector<std::vector<int>> &spatial_neighbor_indices);

 protected:

  int root_rank_;
  int rank_;
  int number_of_mpich_threads_;
  int number_of_elements_per_mpich_thread_;
  std::vector<int> loop_indices_;
  std::vector<int> angular_velocity_loop_indices_;

};

#endif //SPSSFINITEVOLUMEMETHODS_THREAD_HPP
