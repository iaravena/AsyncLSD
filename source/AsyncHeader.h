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

#ifndef ASYNC_INCLUDED
	
	// define ASYNC_INCLUDED to avoid loading the header again
	#define ASYNC_INCLUDED
	
	/* INCLUDE STATEMENTS */
	
	// Standard libraries
	#include <stdio.h>			// provides functions for reading and printing
	#include <stdlib.h>			// provides several functions for managing memory, files and interact with system
	#include <stddef.h>			// definitions (apparently needed to define size_t)
	#include <math.h>			// math functions
	#include <time.h>			// time functions
	
	// MPI
	#include <mpi.h>			// includes MPI library. Use MSMPI in Windows.
	
	// Xpress API
	#include <xprs.h>			// provides functions for managing Xpress from C, C++
	
	// Custom header files
	#include "LowLevelFunc.h"
	#include "Options.h"
	#include "SMPS_Xpress.h"
	
	/* CONSTANTS */
	
	// Decomposition constants
	#define DECOMP_INFINITY				1E+12		// Infty for initiliazing bounds
	#define MAX_WORKERS_PER_SCENARIO	3			// program will stop if there are more than this amounts of workers per scenario
	#define MAX_LENGTH_CRITICAL_QUEUE	5			// maximum lenght of the critical queue, in terms of the number of primal workers
	#define QNORM95						1.644854	// 95% quantile of standard gaussian distribution
	#define NO_UPDATE_SCALING			0.0
	
	// Identifier for non-scenario and non-candidate tasks
	#define NOSCENID					-1
	#define NOCANDID					-2
	#define NO_PREVIOUS_WORKER			-10
	
	// Next job identifiers
	#define NEXT_COMES_DUAL				0
	#define NEXT_COMES_PRIMAL			1
	
	// Candidate status tags
	#define CANDIDATE_FREE_SLOT			0
	#define CANDIDATE_READY				1
	#define CANDIDATE_ACTIVE			2
	#define CANDIDATE_EVALUATED			3
	
	// Task types (used to tag MPI messages)
	#define DUAL_SCEN_INIT				200			// Dual job: (ScenarioNum, x[])
	#define DUAL_SCEN_MILP				201			// Dual job: (ScenarioNum, x[])
	#define DUAL_NONSCEN_MILP			202			// Dual non scenario job: (NOSCENID, LBscen_delayed, msumx_delayed[], msumx_current[], ucenter[])
	#define PRIMAL_FEAS_PROJ			203			// Primal job: (NOCANDID, NOSCENID, u[])
	#define PRIMAL_SCEN_MILP			204			// Primal job: (CandidateID, ScenarioNum, u[])
	
	// Result types (used to tag MPI messages)
	#define RESULT_DUAL_SCEN_INIT		300			// Dual result: (ScenarioNum, XPRS status, worktime, +Infty, LB, u[])
	#define RESULT_DUAL_SCEN_MILP		301			// Dual result: (ScenarioNum, XPRS status, worktime, Obj (UB), LB, u[])
	#define RESULT_DUAL_NONSCEN_MILP	302			// Dual result: (NOSCENID, XPRS status, worktime, Obj (UB), LB stochastic program, u[])
	#define RESULT_PRIMAL_FEAS_PROJ		203			// Primal projection job result: (CandidateID, worktime, u[])
	#define RESULT_PRIMAL_SCEN_MILP		304			// Primal result: (CandidateID, ScenarioNum, XPRS status, worktime, Obj (UB), LB)
	
	// Other messages
	#define RELEASE_WORKER				401
	
	// Worker status
	#define WORKER_COORDINATOR			1000
	#define WORKER_FREE					1001
	#define WORKER_BUSY_DUAL			1002
	#define WORKER_BUSY_PRIMAL			1003
	
	/* STRUCTS */
	typedef struct generic_job
	{
		int task, scenario, candidate;
		void *x_task;
	} generic_job;
	typedef struct slave
	{
		int status;
		double dispatch_timestamp;
		struct generic_job current_job;
	} slave;
	
	/* FUNCTIONS AND PROCEDURES */
	
	// Functions in Coordinator.c
	extern int Coordinator(const int, const int, int*, char***);
	
	// Functions in Workers.c
	extern int ScenarioWorker(const int, int*, char***);
	
	// Functions in HighLevelFunc.c
	extern int read_SMPS(const char*, const char*,
		XPRSprob*, struct string_buffer*, struct string_buffer*,
		struct string_buffer*, int**, int**,
		struct string_buffer*, struct string_buffer*, int**,
		int**, int**, double**, int*, int*,
		int**, double**, double**, const struct options*);
	extern int free_heap_from_read_SMPS(
		XPRSprob*, struct string_buffer*, struct string_buffer*,
		struct string_buffer*, int**, int**,
		struct string_buffer*, struct string_buffer*, int**,
		int**, int**, double**, int**, double**, double**,
		const struct options*);
	extern int load_scenario_subproblem(const char*, const char*, const int,
		const struct string_buffer*, const struct string_buffer*,
		const struct string_buffer*, const int*, const int*,
		const struct string_buffer*, const struct string_buffer*, const int*,
		const int*, const int*, const double*,
		XPRSprob*, struct scen_differences*, const struct options*);
	extern int restore_core_problem(
		const struct string_buffer*, const struct string_buffer*,
		const struct string_buffer*, const int*, const int*,
		const struct string_buffer*, const struct string_buffer*, const int*,
		const int*, const int*, const double*,
		XPRSprob*, struct scen_differences*, const struct options*);
	extern int read_SMPS_coordinator(const char*, const char*,
		struct string_buffer*, double**,
		struct string_buffer*, double **, const struct options*);
	extern int free_heap_from_read_SMPS_coordinator(
		struct string_buffer*, double**,
		struct string_buffer*, double**, const struct options*);
	extern int create_message_types(const int, MPI_Datatype*, MPI_Datatype*,
		MPI_Datatype*, MPI_Datatype*, MPI_Datatype*, MPI_Datatype*);
	extern int free_message_types(MPI_Datatype*, MPI_Datatype*,
		MPI_Datatype*, MPI_Datatype*, MPI_Datatype*, MPI_Datatype*);
	
#endif
