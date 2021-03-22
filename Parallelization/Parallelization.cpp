//
// Created by Nikita Kruk on 24.06.20.
//

#include "Parallelization.hpp"

#if defined(MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY) \
 || defined(MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM) \
 || defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
#endif

#if defined(MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY) \
 || defined(MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM) \
 || defined(MPI_PARAMETER_SCAN)
void LaunchParallelSession(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}
#else
void LaunchParallelSession(int argc, char **argv)
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY) \
 || defined(MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM) \
 || defined(MPI_PARAMETER_SCAN)
void FinalizeParallelSession()
{
  MPI_Finalize();
}
#else
void FinalizeParallelSession()
{

}
#endif