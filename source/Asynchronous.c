/*
Copyright (C) 2020 Ignacio Aravena.

This file is part of AsyncLSD.

AsyncLSD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AsyncLSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AsyncLSD. If not, see <https://www.gnu.org/licenses/>.
*/

/* Asynchronous scenario decomposition with randomized coordinated descent method
 * Version 0.1 (August 1st, 2017)
 * Author: Ignacio Aravena  (ignacio.aravena@uclouvain.be) */

/* Asynchronous algorithm header */
#include "AsyncHeader.h"

int main(int argc, char **argv)
{
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	
	// Check that enough arguments have been passed
	if(argc < 4){
		printf(	"\nAsynchronous scenario decomposition with randomized coordinated"
				"\ndescent method. Version 0.1 (August 1st, 2017)."
				"\n"
				"\nUsage:"
				"\nmpiexec -n $num_procesors AsyncParDecomposition \"$WorkDir\""
				"\n    \"$InstanceDir\" \"$Instance\" \"$OptionsFile\""
				"\n"
				"\n    $num_procesors: Number of MPI processes. This argument might be"
				"\n        passed automatically if your system is managed through Slurm."
				"\n    $WorkDir: Working directoy. Here we will write all files"
				"\n        concerning execution of the algorithm."
				"\n    $InstanceDir: Directory containing the instance in SMPS,"
				"\n        compressed as $Instance.tar.gz."
				"\n    $Instance: Name of the instance to be solved."
				"\n    $OptionsFile: File describing the configuration parameters for"
				"\n        the execution of the decomposition and for the solution of"
				"\n        subproblems using Xpress."
				"\n"
				"\n All arguments must be passed to AsyncParDecomposition.\n");
		return(-1);
	}
	
	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	
	// Print off a hello world message
	printf("\nInitialized process %d/%d - hosted by %s.\n",
		world_rank, world_size, processor_name);
	
	// Check that there are enough workers
	if( world_size < 2 ){
		printf(
			"\n%d processes detected, algorithm requires at least 2 MPI processes"
			"\nto run and it can handle up to %d*(number of scenarios) + 1"
			"\nprocesses. Program will exit.\n", world_size, MAX_WORKERS_PER_SCENARIO);
		return(-1);
	}
	
	// Call routine according to rank
	int outstatus;
	if(world_rank == 0){
		outstatus = Coordinator(world_rank, world_size, &argc, &argv);
		MPI_Abort(MPI_COMM_WORLD, outstatus);
	} else {
		outstatus = ScenarioWorker(world_rank, &argc, &argv);
		if( outstatus != 0 ){
			printf("\nProblem detected in worker %d. Program will exit.\n", world_rank);
			MPI_Abort(MPI_COMM_WORLD, outstatus);
		}
	}
	
	// Finalize the MPI environment
	MPI_Finalize();
	
	// Return success indicator
	return(0);
}