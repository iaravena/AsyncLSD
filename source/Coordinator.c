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

/* 
 * Coordinator routines derivated from main
 */

/* Asynchronous algorithm header */
#include "AsyncHeader.h"
#include "Coordinator.h"

/* Coordinator */
int Coordinator(const int worldrank, const int worldsize, int *argc, char ***argv)
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	// Iterators
	int i, j, k;
	
	// Time registers (BoT: Beggining of Time)
	struct timespec BoT, startt, endt;
	current_utc_time(&BoT);
	
	// Options for execution
	struct options global_opt;
	if( read_options(*((*argv) + 4), &global_opt) != 0 ){
		printf("\nError while reading options. Program will exit.\n");
		return(-1);
	}
	
	// Set seed
	unsigned seed = (unsigned) BoT.tv_nsec;
	srand(seed);
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY)
		printf("\nRandom seed set to %d", seed);
	
	// Prepare files and directory for execution
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY){
		printf("\nCreating execution folder and decompressing SMPS instance ...");
		current_utc_time(&startt);
	}
	if( prepare_files(*((*argv) + 1), *((*argv) + 2), *((*argv) + 3), *((*argv) + 4), &global_opt) != 0 ) {
		printf("\nError while preparing files. Program will exit.\n");
		return(-1);
	}
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY){
		current_utc_time(&endt);
		printf("\nFiles for algorithm execution prepared in %.2f secs.\n",
			((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC);
	}
	
	// Synchronization point: waiting for files to be ready
	MPI_Barrier( MPI_COMM_WORLD );
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY){
		printf("\nStarting Xpress and reading SMPS files...");
		current_utc_time(&startt);
	}
	
	// Variables and pointers for storing problem data
	struct string_buffer first_stage_cols, scenarios;
	double *update_scaling, *probabilities;
	
	// Initialize Xpress-Optimizer
	char *message = (char*) malloc(sizeof(char)*512);
	if( XPRSinit(NULL) ) {
		XPRSgetlicerrmsg(message,512);
		printf("%s\n", message);
	}
	free(message);
	
	// Read first stage variables and scenarios
	if( read_SMPS_coordinator(*((*argv) + 1), *((*argv) + 3),
		&first_stage_cols, &update_scaling,
		&scenarios, &probabilities, &global_opt) != 0 ){
		printf("\nError while reading SMPS files. Program will exit.\n");
		return(-1);
	}
	
	// Terminate Xpress (it is not needed anymore by the Coordinator)
	XPRSfree();
	
	// Synchronization point: waiting for files to be read
	MPI_Barrier( MPI_COMM_WORLD );
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY){
		current_utc_time(&endt);
		printf("\nXpress launched and SMPS files read (by all workers) in %.2f seconds.",
			((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC);
	}
	
	// Check that the number of workers is not more than three times the number of scenarios
	if( worldsize - 1 > MAX_WORKERS_PER_SCENARIO*scenarios.num_elems ){
		printf("\nNumber of workers is more than %d times the number of scenarios."
			"\nProgram will exit in order to avoid inefficient execution.\n", MAX_WORKERS_PER_SCENARIO);
		return(-1);
	}
	// adjust dual share to the number of scenarios
	global_opt.dual_share_start = fmin(global_opt.dual_share_start,
		((double) scenarios.num_elems)/((double) worldsize - 1.0));
	global_opt.dual_share_end = fmin(global_opt.dual_share_end,
		((double) scenarios.num_elems)/((double) worldsize - 1.0));
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	// Define MPI types for communications
	MPI_Datatype dualscenjobmsg, dualnonscenjobmsg, dualjobresult,
		primaljobmsg, primalprojjobresult, primaljobresult;
	if( create_message_types(first_stage_cols.num_elems,
		&dualscenjobmsg, &dualnonscenjobmsg, &dualjobresult,
		&primaljobmsg, &primalprojjobresult, &primaljobresult) != 0 ){
		printf("\nError while creating MPI message types.\n");
		return(-1);
	}
	
	// Allocate space in memory for holding algorithm variables
	int maxnumcandidates = global_opt.max_iterations;
	int primal_queue_size = scenarios.num_elems *
		((int) ceil(2*(1-fmin(global_opt.dual_share_start, global_opt.dual_share_end)) *
			(worldsize-1)/scenarios.num_elems));
	int *candidate_status, *candidate_numevals, *primal_queue_candidate,
		*primal_queue_scenario;
	void *x_evaluated, *v_evaluated, **v_evaluated_scen, *u_evaluated,
		*x_current, **x_current_scen, *x_nonscen, **candidates;
	double *evaluation_timestamp, **x_evaluated_scen_dbl, **v_evaluated_scen_dbl,
		*u_evaluated_dbl, **x_current_scen_dbl, *LBscen, *msumx_evaluated, *msumx_current,
		*ucenter, *candidate_UB, **candidate_UB_scen;
	allocate_memory_space_async_alg(
		scenarios.num_elems, first_stage_cols.num_elems, maxnumcandidates, primal_queue_size,
		&evaluation_timestamp, &x_evaluated, &x_evaluated_scen_dbl, &v_evaluated, &v_evaluated_scen,
		&v_evaluated_scen_dbl, &u_evaluated, &u_evaluated_dbl, &x_current, &x_current_scen,
		&x_current_scen_dbl, &x_nonscen, &LBscen, &msumx_evaluated, &msumx_current, &ucenter,
		&candidates, &candidate_status, &candidate_UB, &candidate_numevals, &candidate_UB_scen,
		&primal_queue_candidate, &primal_queue_scenario);								// here is where all the pointer madness happens
	double LBaux = -DECOMP_INFINITY, LBk = -DECOMP_INFINITY, UBk = DECOMP_INFINITY;		// bound tracking
	double gnorm2;																		// subgradient norm, only computed when evaluating non scenario subproblem
	
	// Initialize dual mutipliers to zero -> OPTION FOR READING MULTIPLIERS TO BE ADDED IN FUTURE VERSIONS
	for(i = 0; i < scenarios.num_elems; i++){
		*((int*) *(x_current_scen + i)) = i;
		for(j = 0; j < first_stage_cols.num_elems; j++){
			*((*(x_current_scen_dbl + i)) + j) = 0.0;
		}
	}
	*((int*) x_nonscen) = NOSCENID;
	*LBscen = -DECOMP_INFINITY;
	for(j = 0; j < first_stage_cols.num_elems; j++){
		*(msumx_current + j) = 0.0;
	}
	
	// Allocate space in memory for parallel execution control variables
	struct slave *workers = (struct slave*) malloc(sizeof(struct slave)*worldsize);
	MPI_Request *work_requests = (MPI_Request*) malloc(sizeof(MPI_Request)*worldsize);
	int numcandidates = 0;
	int primal_queue_numjobs = 0;
	int primal_queue_first, primal_queue_last;
	(*workers).status = WORKER_COORDINATOR;
	for(i = 1; i < worldsize; i++) (*(workers + i)).status = WORKER_FREE;
	
	// Write iterations file header
	FILE *iterfile;
	{
		char *iterfilename = (char*) malloc(sizeof(char)*200);
		strcpy(iterfilename, *((*argv) + 1));
		strcat(iterfilename, SYSTEM_SLASH);
		strcat(iterfilename, "Iterations.csv");
		iterfile = fopen(iterfilename, "w");
		fprintf(iterfile, "Timestamp,Worker,WorkerTime,Iteration,StepType,Scenario,UB,LB,LBBest,DualGap\n");
		free(iterfilename);
	}
	
	// Write candidates file header
	FILE *candidatesfile;
	{
		char *candidatesfilename = (char*) malloc(sizeof(char)*200);
		strcpy(candidatesfilename, *((*argv) + 1));
		strcat(candidatesfilename, SYSTEM_SLASH);
		strcat(candidatesfilename, "Candidates.csv");
		candidatesfile = fopen(candidatesfilename, "w");
		fprintf(candidatesfile, "ID");
		for(i = 0; i < first_stage_cols.num_elems; i++){
			fprintf(candidatesfile, ",%s",
				first_stage_cols.buffer + (first_stage_cols.atom_len+1)*i);
		}
		for(i = 0; i < scenarios.num_elems; i++){
			fprintf(candidatesfile, ",%s",
				scenarios.buffer + (scenarios.atom_len+1)*i);
		}
		fprintf(candidatesfile, ",UB\n");
		free(candidatesfilename);
	}
	
	// Initialize subgradients and lower bounds using period relaxation
	k = -1;
	const int maxCDiterations = scenarios.num_elems*global_opt.max_iterations;
	int num_arrivals, num_projections;
	double av_time_between_arrivals, stddev_time_between_arrivals,
		av_time_projection, stddev_time_projection;
	current_utc_time(&startt);		// here I start counting time of parallel execution
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY)
		printf("\n\nTStmp\tWkr\tWkrT\tIter\tEType\tScen\tUB\tLB\tDualGap\t|g|_2^2\tgTime\n");
	fflush(stdout);
	if( initialize_subgradients_and_bounds(
		&dualscenjobmsg, &dualnonscenjobmsg, &dualjobresult,
		worldsize, workers, work_requests, k, scenarios.num_elems, probabilities,
		first_stage_cols.num_elems, update_scaling, x_current_scen, x_current_scen_dbl,
		x_evaluated_scen_dbl, x_nonscen, LBscen, msumx_evaluated, msumx_current, ucenter,
		evaluation_timestamp, v_evaluated_scen, v_evaluated_scen_dbl,
		u_evaluated, u_evaluated_dbl, &LBaux, &LBk, &gnorm2,
		&numcandidates, &maxnumcandidates, &candidates, &candidate_status, &candidate_UB,
		&candidate_numevals, &candidate_UB_scen,
		&num_arrivals, &av_time_between_arrivals, &stddev_time_between_arrivals,
		&num_projections, &av_time_projection, &stddev_time_projection,
		&startt, iterfile, &global_opt) != 0){
		printf("\nError while initializing algorithm using the period relaxation. Program will exit.\n");
		return(-1);
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	// Launch initial pile of jobs
	int numduals = 0, numprimals = 0, worker_rank = 1;
	const int nonscenario_launching_period =
		(int) ceil(global_opt.nonscen_evaluation_period * scenarios.num_elems);
	int previous_primal_worker = NO_PREVIOUS_WORKER;
	int primal_projection_running = 0;
	int contiguous_scenario_updates = 0;
	int next_candidate;
	k = 0;
	while( worker_rank < worldsize ){
		// launch as many dual jobs as I can
		if( k < scenarios.num_elems ){
			// perform CD update and launch scenario subproblem
			if( update_and_launch_dual_scen_subprob(k, k, gnorm2, LBk, UBk, scenarios.num_elems,
				probabilities, first_stage_cols.num_elems, update_scaling,
				x_current_scen, x_current_scen_dbl,
				msumx_current, v_evaluated_scen_dbl, u_evaluated_dbl, worker_rank, workers,
				work_requests, &dualscenjobmsg, &startt, &global_opt) != 0 ){
				printf("\nError while updating and launching dual scenario subproblem. Program will exit.\n");
				return(-1);
			}
			// increment updates and numduals counters
			k++;
			numduals++;
			contiguous_scenario_updates++;
		// if there are slots available, launch primal jobs
		}else{
			if( global_opt.primal_recovery == PRIMAL_RECOVERY_FIFO
				|| global_opt.primal_recovery == PRIMAL_RECOVERY_RND
				|| global_opt.primal_recovery == PRIMAL_RECOVERY_LIFO ){
				// if there are no pending jobs, activate a new candidate
				if( primal_queue_numjobs == 0 ){
					// look for next primal candidate from the candidate list
					if( select_next_primal_candidate(numcandidates, candidate_status,
						&next_candidate, &global_opt) != 0 ){
						printf("\nError while selecting next candidate from the candidate list. Program will exit.\n");
						return(-1);
					}
					// activate next primal candidate
					if( activate_primal_candidate(next_candidate, scenarios.num_elems,
						numcandidates, candidate_status, candidate_UB_scen,
						&primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
						&primal_queue_last, &primal_queue_candidate, &primal_queue_scenario) != 0 ){
						printf("\nError while activating primal candidate. Program will exit.\n");
						return(-1);
					}
				}
				// launch next primal job from the queue
				if( launch_primal_candidate_evaluation(first_stage_cols.num_elems, numcandidates,
					candidates, primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
					primal_queue_candidate, primal_queue_scenario, worker_rank, &previous_primal_worker,
					workers, work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
					printf("\nError while launching primal candidate evaluation. Program will exit.\n");
					return(-1);
				}
			}else if( global_opt.primal_recovery == PRIMAL_RECOVERY_IS ){
				// launch IS candidate generation and projection
				if( generate_and_project_IS_primal_candidate(
					scenarios.num_elems, probabilities, first_stage_cols.num_elems,
					v_evaluated_scen, v_evaluated_scen_dbl, worker_rank, workers,
					work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
					printf("\nError while generating and launching IS primal candidate. Program will exit.\n");
					return(-1);
				}
				// update primal projection counter
				primal_projection_running++;
			}else{
				printf("\nUnrecognized primal recovery method %d. Program will exit.\n",
					global_opt.primal_recovery);
				return(-1);
			}
			// increment numprimals counter
			numprimals++;
		}
		worker_rank++;
	}
	/*
	2017-04-03: left out to allow serial execution (using 1 Coordinator + 1 Worker)
	// launch non-scenario subproblem and bound computation using the last worker (worker_rank == worldsize - 1)
	if( launch_dual_nonscen_subprob(first_stage_cols.num_elems, x_nonscen,
		worker_rank, workers, work_requests, &dualnonscenjobmsg, &startt,
		&global_opt) != 0 ){
		printf("\nError while launching dual non-scenario subproblems. Program will exit.\n");
		return(-1);
	}
	numduals++;
	*/
	
	// Asynchronous algorithm loop
	double ds = dual_share(k, scenarios.num_elems, &global_opt);
	int primal_queue_critical_lenght = ((int) ceil((1-ds)*(worldsize-1))) - 1;
	double previous_arrival_ts, current_arrival_ts;
	current_utc_time(&endt);
	previous_arrival_ts = ((double) (endt.tv_sec - startt.tv_sec)) +
		((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
	MPI_Status result_status;
	int best_candidate;
	int nonscenario_subproblem_running = 0;
	while( k < maxCDiterations ){
		
		// Probe next message
		if( MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &result_status) != MPI_SUCCESS ){
			printf("\nError while probing for result message (MPI_Probe). Program will exit.\n");
			return(-1);
		}
		
		// Update stats about arrivals
		current_utc_time(&endt);
		current_arrival_ts = ((double) (endt.tv_sec - startt.tv_sec)) +
			((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
		update_time_stats(&num_arrivals, &av_time_between_arrivals, &stddev_time_between_arrivals,
			current_arrival_ts - previous_arrival_ts);
		previous_arrival_ts = current_arrival_ts;
		
		// Dispatch posprocessor depending on the received tag
		switch( result_status.MPI_TAG ){
			fflush(stdout);
			case RESULT_DUAL_SCEN_MILP :
				if( posprocess_dual_scen_subprob(k, scenarios.num_elems, probabilities,
					first_stage_cols.num_elems, evaluation_timestamp, x_evaluated_scen_dbl,
					LBscen, msumx_evaluated, v_evaluated_scen, v_evaluated_scen_dbl, ucenter,
					&result_status, workers, work_requests, &dualjobresult, &numcandidates,
					&maxnumcandidates, &candidates, &candidate_status, &candidate_UB,
					&candidate_numevals, &candidate_UB_scen, &startt, iterfile, &global_opt) != 0 ){
					printf("\nProgram failed to posprocess dual scenario subproblem solution. Program will exit.\n");
					return(-1);
				}
				numduals--;
				break;
			case RESULT_DUAL_NONSCEN_MILP :
			fflush(stdout);
				if( posprocess_dual_nonscen_subprob(k, scenarios.num_elems, probabilities,
					first_stage_cols.num_elems, update_scaling, u_evaluated, u_evaluated_dbl, v_evaluated_scen_dbl,
					UBk, &LBaux, &LBk, &gnorm2, &result_status, workers, work_requests,
					&dualjobresult, &startt, iterfile, &global_opt) != 0 ){
					printf("\nProgram failed to posprocess dual nonscenario subproblem solution. Program will exit.\n");
					return(-1);
				}
				nonscenario_subproblem_running = 0;
				numduals--;
				break;
			case RESULT_PRIMAL_FEAS_PROJ :
				fflush(stdout);
				if( posprocess_primal_feasibility_projection(k, primal_queue_critical_lenght,
					scenarios.num_elems, first_stage_cols.num_elems, &numcandidates, &maxnumcandidates,
					&candidates, &candidate_status, &candidate_UB, &candidate_numevals, &candidate_UB_scen,
					&primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
					&primal_queue_last, &primal_queue_candidate, &primal_queue_scenario,
					&num_projections, &av_time_projection, &stddev_time_projection,
					&result_status, workers, work_requests, &primalprojjobresult,
					&startt, iterfile, &global_opt) != 0){
					printf("\nProgram failed to posprocess feasibility projection subproblem solution. Program will exit.\n");
					return(-1);
				}
				primal_projection_running--;
				numprimals--;
				break;
			case RESULT_PRIMAL_SCEN_MILP :
				if( posprocess_primal_candidate_evaluation(k, scenarios.num_elems, probabilities,
					candidates, candidate_status, candidate_numevals, candidate_UB, candidate_UB_scen,
					&best_candidate, &UBk, LBk, first_stage_cols.num_elems, ucenter, &result_status,
					workers, work_requests, &primaljobresult, &startt, iterfile, candidatesfile,
					&global_opt) != 0 ){
					printf("\nProgram failed to posprocess primal evaluation subproblem solution. Program will exit.\n");
					return(-1);
				}
				numprimals--;
				break;
			default :
				printf("\nUnrecognized tag %d received. Program will exit.\n", result_status.MPI_TAG);
				return(-1);
		}
		fflush(stdout);
		
		// Check tolerance
		if( UBk - LBk < 0 ){
			printf("\nUpper bound %f is smaller than lower bound %f. Program will exit.\n", UBk, LBk);
			return(-1);
		}else if( UBk - LBk <= global_opt.rel_dual_gap_tol*(1E-10+fmax(fabs(UBk),fabs(LBk))) ){
			printf("\nTarget relative duality gap achieved. Exiting main loop.\n");
			break;
		}
		
		// Check max time
		current_utc_time(&endt);
		if( ((double) (endt.tv_sec - BoT.tv_sec)) +
				((double) (endt.tv_nsec - BoT.tv_nsec))*NANOSEC2SEC >= global_opt.max_time ){
			printf("\nAlgorithm hitted the time limit. Exiting main loop.\n");
			break;
		}
		
		// Decide what to do next based on the current assignement of workers and the dual share
		//printf("\nI'm here, next job will be %d (%d, %d, %d, %d, %f).\n",
		//	next_job(scenarios.num_elems, worldsize - 1, numduals, numprimals, ds),
		//	scenarios.num_elems, worldsize - 1, numduals, numprimals, ds);
		worker_rank = result_status.MPI_SOURCE;
		if( next_job(scenarios.num_elems, worldsize - 1, numduals, numprimals, ds) == NEXT_COMES_DUAL ){
			// Launch dual job
			if( contiguous_scenario_updates < nonscenario_launching_period
				|| nonscenario_subproblem_running == 1){
				// perform CD update and launch scenario subproblem
				if( update_and_launch_dual_scen_subprob(k, RANDOM_SCENARIO, gnorm2, LBaux, UBk, scenarios.num_elems,
					probabilities, first_stage_cols.num_elems, update_scaling, x_current_scen, x_current_scen_dbl,
					msumx_current, v_evaluated_scen_dbl, u_evaluated_dbl, worker_rank, workers,
					work_requests, &dualscenjobmsg, &startt, &global_opt) != 0 ){
					printf("\nError while updating and launching dual scenario subproblem. Program will exit.\n");
					return(-1);
				}
				// increment updates counters
				k++;
				contiguous_scenario_updates++;
				// write current multipliers every S updates
				if(k % scenarios.num_elems == 0){
					write_dual_multipliers(*((*argv) + 1), &scenarios, &first_stage_cols, x_current_scen_dbl);
				}
				// recompute dual share
				ds = dual_share(k, scenarios.num_elems, &global_opt);
			}else{
				// launch non-scenario subproblem and bound computation
				if( launch_dual_nonscen_subprob(first_stage_cols.num_elems, x_nonscen,
					worker_rank, workers, work_requests, &dualnonscenjobmsg, &startt,
					&global_opt) != 0 ){
					printf("\nError while launching dual non-scenario subproblems. Program will exit.\n");
					return(-1);
				}
				// update non scenario subproblem execution status
				contiguous_scenario_updates = 0;
				nonscenario_subproblem_running = 1;
			}
			numduals++;
		}else{
			// Launch primal job
			if( global_opt.primal_recovery == PRIMAL_RECOVERY_FIFO
				|| global_opt.primal_recovery == PRIMAL_RECOVERY_RND
				|| global_opt.primal_recovery == PRIMAL_RECOVERY_LIFO ){
				// compute critical queue length
				primal_queue_critical_lenght = (int) ceil((1-ds)*(worldsize-1));
				if( primal_queue_numjobs <= primal_queue_critical_lenght ){
					// look for next primal candidate from the candidate list
					if( select_next_primal_candidate(numcandidates, candidate_status,
						&next_candidate, &global_opt) != 0 ){
						printf("\nError while selecting next candidate from the candidate list. Program will exit.\n");
						return(-1);
					}
					// activate next primal candidate
					if( activate_primal_candidate(next_candidate, scenarios.num_elems,
						numcandidates, candidate_status, candidate_UB_scen,
						&primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
						&primal_queue_last, &primal_queue_candidate, &primal_queue_scenario) != 0 ){
						printf("\nError while activating primal candidate. Program will exit.\n");
						return(-1);
					}
				}
				// launch next primal job from the queue
				if( launch_primal_candidate_evaluation(first_stage_cols.num_elems, numcandidates,
					candidates, primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
					primal_queue_candidate, primal_queue_scenario, worker_rank, &previous_primal_worker,
					workers, work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
					printf("\nError while launching primal candidate evaluation. Program will exit.\n");
					return(-1);
				}
			}else if( global_opt.primal_recovery == PRIMAL_RECOVERY_IS ){
				primal_queue_critical_lenght = (int) fmax( ceil((1-ds)*(worldsize-1)),
					time_critical_queue_length(av_time_projection, stddev_time_projection,
						num_arrivals, av_time_between_arrivals, stddev_time_between_arrivals,
						(int) ceil((1-ds)*(worldsize-1))));
				if( primal_queue_numjobs > primal_queue_critical_lenght
					|| (primal_queue_numjobs <= primal_queue_critical_lenght
						&& primal_queue_numjobs > 0
						&& primal_projection_running > 0) ){
					// launch next primal job from the queue
					if( launch_primal_candidate_evaluation(first_stage_cols.num_elems, numcandidates,
						candidates, primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
						primal_queue_candidate, primal_queue_scenario, worker_rank, &previous_primal_worker,
						workers, work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
						printf("\nError while launching primal candidate evaluation. Program will exit.\n");
						return(-1);
					}
				}else if(primal_queue_numjobs <= primal_queue_critical_lenght
					&& primal_queue_numjobs > 0
					&& primal_projection_running == 0){
					// launch IS candidate generation and projection
					if( generate_and_project_IS_primal_candidate(
						scenarios.num_elems, probabilities, first_stage_cols.num_elems,
						v_evaluated_scen, v_evaluated_scen_dbl, worker_rank, workers,
						work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
						printf("\nError while generating and launching IS primal candidate. Program will exit.\n");
						return(-1);
					}
					// update primal projection counter
					primal_projection_running++;
				}else if( primal_queue_numjobs == 0 ){
					// if there are no jobs left, try to activate the last one ready in the candidates list
					for(i = numcandidates - 1; i >= 0; i--){
						if( *(candidate_status + i) == CANDIDATE_READY ) break;
					}
					if(i >= 0){ // there is a candidate -> activate it and launch primal evaluation 
						if( activate_primal_candidate(i, scenarios.num_elems,
							numcandidates, candidate_status, candidate_UB_scen,
							&primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
							&primal_queue_last, &primal_queue_candidate, &primal_queue_scenario) != 0 ){
							printf("\nError while activating primal candidate. Program will exit.\n");
							return(-1);
						}
						if( launch_primal_candidate_evaluation(first_stage_cols.num_elems, numcandidates,
							candidates, primal_queue_size, &primal_queue_numjobs, &primal_queue_first,
							primal_queue_candidate, primal_queue_scenario, worker_rank, &previous_primal_worker,
							workers, work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
							printf("\nError while launching primal candidate evaluation. Program will exit.\n");
							return(-1);
						}
					}else{ // there is no READY candidate -> launch IS candidate generation and projection
						if( generate_and_project_IS_primal_candidate(
							scenarios.num_elems, probabilities, first_stage_cols.num_elems,
							v_evaluated_scen, v_evaluated_scen_dbl, worker_rank, workers,
							work_requests, &primaljobmsg, &startt, &global_opt) != 0 ){
							printf("\nError while generating and launching IS primal candidate. Program will exit.\n");
							return(-1);
						}
						// update primal projection counter
						primal_projection_running++;					}
				}else{
					printf("\nUnrecognized combination of algorithm parameters for primal execution."
						"\nProgram will exit.\n");
					return(-1);
				}
			}
			numprimals++;
		}
	}
	
	// Release workers
	MPI_Request dummy_req;
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY)
		printf("\n\nSending termination signal to all workers...");
	for(i = 0; i < worldsize; i++){
		if( i != worldrank ){
			MPI_Isend(NULL, 0, MPI_INT, i, RELEASE_WORKER, MPI_COMM_WORLD, &dummy_req);
			MPI_Request_free(&dummy_req);
		}
	}
	
	// Write final multipliers
	write_dual_multipliers(*((*argv) + 1), &scenarios, &first_stage_cols, x_current_scen_dbl);
	
	// Close output to files
	fclose(iterfile);
	fclose(candidatesfile);
	
	// Free heap
	free_memory_candidates(numcandidates, candidate_status, candidates, candidate_UB_scen);
	free_message_types(&dualscenjobmsg, &dualnonscenjobmsg, &dualjobresult,
		&primaljobmsg, &primalprojjobresult, &primaljobresult);
	free_memory_space_async_alg(&evaluation_timestamp, &x_evaluated, &x_evaluated_scen_dbl,
		&v_evaluated, &v_evaluated_scen, &v_evaluated_scen_dbl, &u_evaluated, &x_current,
		&x_current_scen, &x_current_scen_dbl, &x_nonscen,
		&candidates, &candidate_status, &candidate_UB, &candidate_numevals, &candidate_UB_scen,
		&primal_queue_candidate, &primal_queue_scenario);
	if( free_heap_from_read_SMPS_coordinator(
		&first_stage_cols, &update_scaling,
		&scenarios, &probabilities, &global_opt) != 0 ){
		printf("\nError while cleaning the heap from variables storing the SMPS instance. Program will exit.\n");
		return(-1);
	}
	free_options(&global_opt);
	
	// Final message
	current_utc_time(&endt);
	if(global_opt.verbosity >= VERB_COORDINATOR_ONLY){
		printf("\n\nFinal UB:\t%f\t(incumbent candidate %d)"
			"\nFinal LB:\t%f"
			"\nFinal gap:\t%f%%\n",
			UBk, best_candidate, LBk, 100*(UBk - LBk)/(1E-10+fmax(fabs(UBk),fabs(LBk))) );
		printf("\n\nAll done. Elapsed time:\t%.2f secs\n\n",
			((double) (endt.tv_sec - BoT.tv_sec)) +
				((double) (endt.tv_nsec - BoT.tv_nsec))*NANOSEC2SEC);
	}
	
	// Return success indicator
	return(0);
}
