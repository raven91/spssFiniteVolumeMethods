//
// Created by Nikita Kruk on 24.06.20.
//

#ifndef SPSSFINITEVOLUMEMETHODS_PARALLELIZATION_HPP
#define SPSSFINITEVOLUMEMETHODS_PARALLELIZATION_HPP

#include "../Definitions.hpp"

#if defined(MPI_PARALLELIZATION_IN_OMEGA_ONE_SIDED_SHARED_MEMORY) \
 || defined(MPI_PARALLELIZATION_IN_ALL_ONE_SIDED_SHARED_MEMORY_RANDOM) \
 || defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
const MPI_Datatype kRealTypeForMpi = MPI_DOUBLE;
#else
const int kRealTypeForMpi = -1;
#endif

void LaunchParallelSession(int argc, char **argv);
void FinalizeParallelSession();

#endif //SPSSFINITEVOLUMEMETHODS_PARALLELIZATION_HPP
