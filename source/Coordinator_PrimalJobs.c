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
 * Routines derived from Coordinator
 */
#include "AsyncHeader.h"
#include "Coordinator.h"

//	*primal_proj_receiver = malloc( sizeof(int) + sizeof(double)*(1+numfirststagecols) );			// (candj, worker time, feasible candidate[])
//	*primal_cand_scen_receiver = malloc( sizeof(int)*3 + sizeof(double)*3 );						// (candj, sceni, XPRS status, worker time, fi(candj)_UB, fi(candj)_LB)

/* Launch next primal evaluation job in the queue */
int launch_primal_candidate_evaluation(
	const int numfirststagecols, const int numcandidates, void **candidates,
	const int primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
	const int *primal_queue_candidate, const int *primal_queue_scenario,
	const int assigned_worker, int *previous_primal_worker,
	struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *primaljobmsg, struct timespec *startt, const struct options *opt)
{
	// Check whether the assigned worker is actually free
	if( (*(workers + assigned_worker)).status != WORKER_FREE ){
		printf("\nWorker assigned for next primal candidate evaluation is not free."
			"\nWorker status: %d (it should be %d)\n",
			(*(workers + assigned_worker)).status, WORKER_FREE);
		return(-1);
	}
	
	// Check that there are jobs left
	if( *primal_queue_numjobs == 0 ){
		printf("\nNo jobs left in the primal queue.\n");
		return(-1);
	}
	
	// Wait for previous primal MPI_Isend to complete
	// this is done in order to avoid modifying the scenario dummy in the candidates buffer
	if( *previous_primal_worker != NO_PREVIOUS_WORKER ){
		MPI_Status send_status;
		MPI_Wait(work_requests + *previous_primal_worker, &send_status);
	}
	
	// Allocate task details in the workers object
	struct timespec endt;
	current_utc_time(&endt);
	(*(workers + assigned_worker)).status = WORKER_BUSY_PRIMAL;
	(*(workers + assigned_worker)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
	(*(workers + assigned_worker)).current_job.task = PRIMAL_SCEN_MILP;
	(*(workers + assigned_worker)).current_job.candidate = *(primal_queue_candidate + *primal_queue_first);
	(*(workers + assigned_worker)).current_job.scenario = *(primal_queue_scenario + *primal_queue_first);
	*(((int*) *(candidates + (*(workers + assigned_worker)).current_job.candidate)) + 1) =
		(*(workers + assigned_worker)).current_job.scenario;
	
	// Launch primal job
	if( MPI_Isend(*(candidates + (*(workers + assigned_worker)).current_job.candidate), 1, *primaljobmsg,
		assigned_worker, (*(workers + assigned_worker)).current_job.task, MPI_COMM_WORLD,
		work_requests + assigned_worker) != MPI_SUCCESS ){
		printf("\nError while sending primal evalation subproblem job through MPI_Isend.\n");
		return(-1);
	}
	
	// Update first position and number of primal jobs in the queue
	*primal_queue_first = cycling_next_element(primal_queue_size, *primal_queue_first);
	*primal_queue_numjobs -= 1;
	
	// Update previous primal worker
	*previous_primal_worker = assigned_worker;
	
	// Return success indicator
	return(0);
}

/* Posprocessing a completed primal candidate evaluation for a single scenario */
int posprocess_primal_candidate_evaluation(const int k,
	const int numscenarios, const double *scenprobabilities,
	void **candidates, int *candidate_status, int *candidate_numevals,
	double *candidate_UB, double **candidate_UB_scen,
	int *best_candidate, double *UBk, const double LBk,
	const int numfirststagecols, double *ucenter,
	MPI_Status *result_status, struct slave *workers,
	MPI_Request *work_requests, MPI_Datatype *primaljobresult,
	struct timespec *startt, FILE *iterfile, FILE *candidatefile, const struct options *opt)
{
	// Define internal variables
	int incumbent_worker = (*result_status).MPI_SOURCE;
	int incumbent_candidate = (*(workers + incumbent_worker)).current_job.candidate;
	int incumbent_scenario = (*(workers + incumbent_worker)).current_job.scenario;
	
	// Receive the message
	// Primal result: (CandidateID, ScenarioNum, XPRS status, worktime, Obj, LB)
	void *primal_cand_scen_receiver = malloc( sizeof(int)*3 + sizeof(double)*3 );
	if( MPI_Recv(primal_cand_scen_receiver, 1, *primaljobresult,
		(*result_status).MPI_SOURCE, (*result_status).MPI_TAG,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
		printf("\nError while receiving message %d from worker %d.\n",
			(*result_status).MPI_TAG, (*result_status).MPI_SOURCE);
		return(-1);
	}
	
	// Check that the received scenario coincides with the one registered by the worker
	if( incumbent_candidate != *((int*) primal_cand_scen_receiver)
		|| incumbent_scenario != *(((int*) primal_cand_scen_receiver) + 1) ){
		printf("\nInconsistent reception of primal evaluation."
			"\nExpected (candidate,scenario) (%d,%d), received (%d,%d).\n",
			incumbent_candidate, incumbent_scenario, *((int*) primal_cand_scen_receiver),
			*(((int*) primal_cand_scen_receiver) + 1));
		return(-1);
	}
	
	// Print received message
	struct timespec endt;
	current_utc_time(&endt);
	if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
		printf("%5.1f\t%d\t%5.1f\t%d\tPRIM\t(%d,%d)\t%5.2f\t%5.2f\t-----\t-----\t-----\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			incumbent_worker,
			*((double*) (((int*) primal_cand_scen_receiver) + 3)),
			k,
			*((int*) primal_cand_scen_receiver), *(((int*) primal_cand_scen_receiver) + 1),
			*(((double*) (((int*) primal_cand_scen_receiver) + 3)) + 1)/1E6,
			*(((double*) (((int*) primal_cand_scen_receiver) + 3)) + 2)/1E6);
	// write iteration data
	fprintf(iterfile, "%f,%d,%f,%d,PRIMAL,\"(%d,%d)\",%f,%f,,\n",
		((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
		incumbent_worker,
		*((double*) (((int*) primal_cand_scen_receiver) + 3)),
		k,
		*((int*) primal_cand_scen_receiver), *(((int*) primal_cand_scen_receiver) + 1),
		*(((double*) (((int*) primal_cand_scen_receiver) + 3)) + 1)/1E6,
		*(((double*) (((int*) primal_cand_scen_receiver) + 3)) + 2)/1E6);
	
	// Place the upper bound in the corresponding buffer and increase the count of evaluations
	*((*(candidate_UB_scen + incumbent_candidate)) + incumbent_scenario) = 
		*(((double*) (((int*) primal_cand_scen_receiver) + 3)) + 1);
	*(candidate_numevals + incumbent_candidate) += 1;
	
	// Check whether we completed the evaluation of this candidate
	if( *(candidate_numevals + incumbent_candidate) == numscenarios ){
		// Candidate is evaluated for all scenarios -> consolidate
		consolidate_primal_candidate(incumbent_candidate, numscenarios,
			scenprobabilities, numfirststagecols, candidates, candidate_status,
			candidate_UB, candidate_UB_scen, candidatefile);
		// Check whether we have found a better primal solution
		if( *(candidate_UB + incumbent_candidate) < *UBk ){
			// Update upper bound for the stochastic problem
			*best_candidate = incumbent_candidate;
			*UBk = *(candidate_UB + incumbent_candidate);
			// Update center for smoothing
			if( (*opt).nonscen_center_type == CENTER_INCUMBENT ){
				memcpy((void*) ucenter, (void*) (((int*) *(candidates + incumbent_candidate)) + 2),
					sizeof(double)*numfirststagecols);
			}
			// Print new upper bound message
			double dgap = 100*(*UBk - LBk)/(1E-10 + fmax(fabs(*UBk),fabs(LBk)));
			current_utc_time(&endt);
			if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
				printf("%5.1f\t-----\t-----\t%d\tPRIM\t(%d,All)\t%5.2f\t%5.2f\t%5.2f%%\t-----\t-----\n",
					((double) (endt.tv_sec - (*startt).tv_sec)) +
						((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
					k, incumbent_candidate, *UBk/1E6, LBk/1E6, dgap);
			// write iteration data
			fprintf(iterfile, "%f,,,%d,PRIMAL,\"(%d,All)\",%f,%f,,%f\n",
				((double) (endt.tv_sec - (*startt).tv_sec)) +
					((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
				k, incumbent_candidate, *UBk/1E6, LBk/1E6, dgap);
		}
	}
	
	// Free receiver buffer space
	free(primal_cand_scen_receiver);
	
	// Release worker
	(*(workers + incumbent_worker)).status = WORKER_FREE;
	
	// Return success indicator
	return(0);
}

/* Generate a new IS candidate and launch feasibility projection */
int generate_and_project_IS_primal_candidate(
	const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, void **v_evaluated_scen, double **v_evaluated_scen_dbl,	
	const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *primaljobmsg, struct timespec *startt, const struct options *opt)
{
	// Check whether the assigned worker is actually free
	if( (*(workers + assigned_worker)).status != WORKER_FREE ){
		printf("\nWorker assigned for IS feasibility projection is not free."
			"\nWorker status: %d (it should be %d)\n",
			(*(workers + assigned_worker)).status, WORKER_FREE);
		return(-1);
	}
	
	// Generate IS infeasible candidate
	(*(workers + assigned_worker)).current_job.x_task =
		malloc(sizeof(int)*2 + sizeof(double)*numfirststagecols);
	double *new_candidate =  (double*) (((int*) (*(workers + assigned_worker)).current_job.x_task) + 2);
	generate_new_IS_primal_candidate(numscenarios, scenprobabilities, numfirststagecols,
		v_evaluated_scen, v_evaluated_scen_dbl, new_candidate, opt);
	
	// Allocate task details in the workers object
	struct timespec endt;
	current_utc_time(&endt);
	(*(workers + assigned_worker)).status = WORKER_BUSY_PRIMAL;
	(*(workers + assigned_worker)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
	(*(workers + assigned_worker)).current_job.task = PRIMAL_FEAS_PROJ;
	(*(workers + assigned_worker)).current_job.candidate = NOCANDID;
	(*(workers + assigned_worker)).current_job.scenario = NOSCENID;
	*((int*) (*(workers + assigned_worker)).current_job.x_task) = NOCANDID;
	*(((int*) (*(workers + assigned_worker)).current_job.x_task) + 1) = NOSCENID;
	
	// Launch projection job
	if( MPI_Isend((*(workers + assigned_worker)).current_job.x_task, 1, *primaljobmsg,
		assigned_worker, (*(workers + assigned_worker)).current_job.task, MPI_COMM_WORLD,
		work_requests + assigned_worker) != MPI_SUCCESS ){
		printf("\nError while sending feasibility projection subproblem job through MPI_Isend.\n");
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Posprocessing a completed primal feasibility projection */
int posprocess_primal_feasibility_projection(const int k, const int critical_queue_lenght,
	const int numscenarios, const int numfirststagecols,
	int *numcandidates, int *maxnumcandidates, void ***candidates, int **candidate_status,
	double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
	int *primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
	int *primal_queue_last, int **primal_queue_candidate, int **primal_queue_scenario,
	int *num_projections, double *av_time_projection, double *stddev_time_projection,
	MPI_Status *result_status, struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *primalprojjobresult, struct timespec *startt, FILE *iterfile,
	const struct options *opt)
{
	// Get worker number
	int incumbent_worker = (*result_status).MPI_SOURCE;
	void *primal_proj_receiver = malloc( sizeof(int) + sizeof(double)*(1+numfirststagecols) );
	double *new_candidate_feasible = ((double*) (((int*) primal_proj_receiver) + 1)) + 1;
	
	// Receive message
	if( MPI_Recv(primal_proj_receiver, 1, *primalprojjobresult,
		incumbent_worker, RESULT_PRIMAL_FEAS_PROJ,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
		printf("\nError while receiving message %d from worker %d.\n",
			RESULT_PRIMAL_FEAS_PROJ, incumbent_worker);
		return(-1);
	}
	
	// Check that received message corresponds to projection
	if( *((int*) primal_proj_receiver) != NOCANDID ){
		printf("\nInconsistent reception of subproblem solution. Expected no-candidate %d, received %d.\n",
			NOCANDID, *((int*) primal_proj_receiver));
		return(-1);
	}
	
	// Print received message
	struct timespec endt;
	current_utc_time(&endt);
	if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
		printf("%5.1f\t%d\t%5.1f\t%d\tPROJ\t-----\t-----\t-----\t-----\t-----\t-----\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			incumbent_worker, *((double*) (((int*) primal_proj_receiver) + 1)), k);
	
	// Write iteration data
	fprintf(iterfile, "%f,%d,%f,%d,PROJ,,,,,\n",
		((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
		incumbent_worker, *((double*) (((int*) primal_proj_receiver) + 1)), k);
	
	// Update projection statistics
	update_time_stats(num_projections, av_time_projection, stddev_time_projection,
		*((double*) (((int*) primal_proj_receiver) + 1)));
	
	// Add candidate to the list and activate it
	if( (*opt).primal_recovery == PRIMAL_RECOVERY_IS ){
		int add_status = add_primal_candidate(
			numscenarios, numfirststagecols, new_candidate_feasible,
			numcandidates, maxnumcandidates, candidates, candidate_status,
			candidate_UB, candidate_numevals, candidate_UB_scen, opt);
		if(add_status < 0){
			printf("\nError adding a new primal candidate into the list.\n");
			return(-1);
		}
		// activate candidate only if:
		if( add_status == 0										// the candidate is new
			&& *primal_queue_numjobs <= critical_queue_lenght	// we are at or below the critical primal queue length
			&& activate_primal_candidate(*numcandidates - 1, numscenarios, *numcandidates,
				*candidate_status, *candidate_UB_scen, primal_queue_size, primal_queue_numjobs,
				primal_queue_first, primal_queue_last, primal_queue_candidate,
				primal_queue_scenario) != 0){
			printf("\nError activating primal candidate.\n");
			return(-1);
		}
	}else{
		printf("\nUsing feasibility projection for non IS primal recovery %d.\n",
			(*opt).primal_recovery);
		return(-1);
	}
	
	// Release worker
	(*(workers + incumbent_worker)).status = WORKER_FREE;
	free((*(workers + incumbent_worker)).current_job.x_task);
	
	// Free primal receiver
	free(primal_proj_receiver);
	
	// Return success indicator
	return(0);
}

/* Select next primal candidate from the list */
int select_next_primal_candidate(const int numcandidates, int *candidate_status,
	int *next_candidate_ID, const struct options *opt)
{
	// Select next candidate based on the specified options
	int i, j, *ready_candidates;
	switch( (*opt).primal_recovery ){
		case PRIMAL_RECOVERY_FIFO :
			// look for the first READY candidate
			for(i = 0; i < numcandidates; i++){
				if( *(candidate_status + i) == CANDIDATE_READY ){
					break;
				}
			}
			if( i < numcandidates ){
				*next_candidate_ID = i;
			}else{
				printf("\nThere are no new primal candidates to explore in the candidate list.\n");
				return(-1);
			}
			break;
		case PRIMAL_RECOVERY_RND :
			// enumerate READY candidates, then pick one at random
			j = 0;
			ready_candidates = (int*) malloc(sizeof(int)*numcandidates);
			for(i = 0; i < numcandidates; i++){
				if( *(candidate_status + i) == CANDIDATE_READY ){
					*(ready_candidates + j) = i;
					j++;
				}
			}
			if( j > 0 ){
				*next_candidate_ID = *(ready_candidates + (rand() % j));
			}else{
				printf("\nThere are no new primal candidates to explore in the candidate list.\n");
				return(-1);
			}
			free(ready_candidates);
			break;
		case PRIMAL_RECOVERY_LIFO :
			// look for the last READY candidate
			for(i = numcandidates-1; i >= 0; i--){
				if( *(candidate_status + i) == CANDIDATE_READY ){
					break;
				}
			}
			if( i >= 0 ){
				*next_candidate_ID = i;
			}else{
				printf("\nThere are no new primal candidates to explore in the candidate list.\n");
				return(-1);
			}
			break;
		default:
			printf("\nUnrecognized primal recovery method %d.\n", (*opt).primal_recovery);
			return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Create new primal candidate using the Importance Sampling heuristic */
void generate_new_IS_primal_candidate(const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
	double *new_candidate, const struct options *opt)
{
	// Iterators
	int i, j;
	
	// Compute importances
	double *importance = (double*) malloc(sizeof(double)*numscenarios);
	for(i = 0; i < numscenarios; i++){
		// (sceni, XPRS status, worker time, fi_UB, fi_LB, vevali[]
		*(importance + i) = *(scenprobabilities + i) *
			*(((double*) (((int*) *(v_evaluated_scen + i)) + 2)) + 2);
	}
	
	// Get a sample based on the importances
	int samplesize = (int) ceil((*opt).primal_IS_ssize*numscenarios);
	int tablesize = 0;
	int *samplescenarios = (int*) malloc(sizeof(int)*numscenarios);
	int *frequencies = (int*) malloc(sizeof(int)*numscenarios);
	nonuniformsample(numscenarios, importance, samplesize,
		samplescenarios, frequencies, &tablesize);
	
	// Average the samples to get the candidate
	double weight;
	for(i = 0; i < numfirststagecols; i++){
		*(new_candidate + i) = 0.0;
	}
	for(i = 0; i < tablesize; i++){
		weight = ((double) *(frequencies + i))/((double) samplesize);
		for(j = 0; j < numfirststagecols; j++){
			*(new_candidate + j) += weight * *((*(v_evaluated_scen_dbl + *(samplescenarios + i))) + j);
		}
	}
	
	// Release allocated memory
	free(importance);
	free(samplescenarios);
	free(frequencies);
}

/* Add new primal candidate to the queue. new_candidate must be feasible.
 * Returns: 0 if new_candidate was added, 1 if candidate already existed and -1 if there is an error */
int add_primal_candidate(const int numscenarios,
	const int numfirststagecols, const double *new_candidate,
	int *numcandidates, int *maxnumcandidates, void ***candidates, int **candidate_status,
	double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
	const struct options *opt)
{
	// Check that new_candidate is not already among candidates (backwards)
	{
		int i;
		double *candidate_dbl;
		for(i = (*numcandidates) - 1; i >= 0 ; i--){
			// *candidates {(candj, scenario dummy, uj[]) : j=1,...,numcandidates}
			candidate_dbl = (double*) (((int*) *((*candidates) + i)) + 2);
			if( L1norm(numfirststagecols, new_candidate, candidate_dbl) <  DBL_EPSILON ){
				return(1);	// candidate already considered
			}
		}
	}
	
	// Realloc, if I ran out of space
	if( *numcandidates == *maxnumcandidates ){
		int i, newcandlistsize = (*maxnumcandidates) + numscenarios;
		*candidates = (void**) realloc( *candidates, sizeof(void*) * newcandlistsize );
		*candidate_status = (int*) realloc( *candidate_status, sizeof(int) * newcandlistsize );			// {statusj : j=1,...,numcandidates}
		*candidate_UB = (double*) realloc( *candidate_UB, sizeof(double) * newcandlistsize );			// {UBj : j=1,...,numcandidates}
		*candidate_numevals = (int*) realloc( *candidate_numevals, sizeof(int) * newcandlistsize );		// {num evaluated scenarios j : j=1,...,numcandidates}
		*candidate_UB_scen = (double**) realloc( *candidate_UB_scen, sizeof(double*) * newcandlistsize );
		if( *candidates == NULL || *candidate_status == NULL || *candidate_UB == NULL ||
			*candidate_numevals == NULL || *candidate_UB_scen == NULL ){
			printf("\nrealloc failed to allocate more memory for holding primal candidates.\n");
			return(-1);
		}
		for(i = *maxnumcandidates; i < newcandlistsize; i++){
			*((*candidate_status) + i) = CANDIDATE_FREE_SLOT;
		}
		*maxnumcandidates = newcandlistsize;
	}
	
	// Add new_candidate : (candj, scenario dummy, uj[])
	*((*candidates) + *numcandidates) = malloc( sizeof(int)*2 + sizeof(double)*numfirststagecols );
	*((int*) *((*candidates) + *numcandidates)) = *numcandidates;
	*(((int*) *((*candidates) + *numcandidates)) + 1) = NOSCENID;
	memcpy((void*) (((int*) *((*candidates) + *numcandidates)) + 2), (void*) new_candidate,
		sizeof(double)*numfirststagecols);
	*((*candidate_status) + *numcandidates) = CANDIDATE_READY;
	*((*candidate_UB) + *numcandidates) = DECOMP_INFINITY;
	*((*candidate_numevals) + *numcandidates) = 0;
	*numcandidates += 1;
	
	// Return success indicator
	return(0);
}

/* Activate primal candidate */
int activate_primal_candidate(const int candidate_ID, const int numscenarios,
	const int numcandidates, int *candidate_status, double **candidate_UB_scen,
	int *primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
	int *primal_queue_last, int **primal_queue_candidate, int **primal_queue_scenario)
{
	// Iterators
	int i, j;
	
	// Check that the candidate is ready
	if( candidate_ID >= numcandidates || *(candidate_status + candidate_ID) != CANDIDATE_READY ){
		printf("\nCandidate %d has status %d, it should have status %d to be activated.\n",
			candidate_ID, *(candidate_status + candidate_ID), CANDIDATE_READY);
		return(-1);
	}
	
	// Realloc if I ran out of space in the primal queue
	if( (*primal_queue_size) - (*primal_queue_numjobs) < numscenarios ){
		// create new, larger, queue
		int new_queue_size = *primal_queue_size + numscenarios;
		int *new_queue_candidate = (int*) malloc(sizeof(int)*new_queue_size);
		int *new_queue_scenario = (int*) malloc(sizeof(int)*new_queue_size);
		// assign elements of the current queue into the new queue
		i = 0;
		j = *primal_queue_first;
		for(i = 0; i < *primal_queue_numjobs; i++){
			*(new_queue_candidate + i) = *((*primal_queue_candidate) + j);
			*(new_queue_scenario + i) = *((*primal_queue_scenario) + j);
			j = cycling_next_element(*primal_queue_size, j);
		}
		// free current memory
		free(*primal_queue_candidate);
		free(*primal_queue_scenario);
		// assign the new memory space to the queue variables (numjobs remain constant)
		*primal_queue_candidate = new_queue_candidate;
		*primal_queue_scenario = new_queue_scenario;
		*primal_queue_first = 0;
		*primal_queue_last = (*primal_queue_numjobs) - 1;
		*primal_queue_size = new_queue_size;
	}
	
	// Add jobs to the queue
	if( *primal_queue_numjobs == 0 ){
		*primal_queue_first = 0;	// start filling from the top
		j = (*primal_queue_size) - 1;
	}else{
		j = *primal_queue_last;
	}
	for(i = 0; i < numscenarios; i++){
		j = cycling_next_element(*primal_queue_size, j);
		*((*primal_queue_candidate) + j) = candidate_ID;
		*((*primal_queue_scenario) + j) = i;
	}
	*primal_queue_last = j;
	*primal_queue_numjobs += numscenarios;
	
	// Allocate space for holding scenario upper bounds
	*(candidate_UB_scen + candidate_ID) = (double*) malloc(sizeof(double)*numscenarios);
	
	// Change candidate status
	*(candidate_status + candidate_ID) = CANDIDATE_ACTIVE;
	
	// Return success indicator
	return(0);
}

/* Consolidate primal candidate solution */
void consolidate_primal_candidate(const int candidate_ID,
	const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, void **candidates, int *candidate_status,
	double *candidate_UB, double **candidate_UB_scen, FILE *candidatesfile)
{
	// Compute candidate's upper bound
	int i;
	*(candidate_UB + candidate_ID) = 0.0;
	for(i = 0; i < numscenarios; i++){
		*(candidate_UB + candidate_ID) += *(scenprobabilities + i) *
			*((*(candidate_UB_scen + candidate_ID)) + i);
	}
	
	// Write candidate to file, along with the bounds for each scenario
	double *candidate_dbl = (double*) (((int*) *(candidates + candidate_ID)) + 2);
	fprintf(candidatesfile, "%d", candidate_ID);
	for(i = 0; i < numfirststagecols; i++){
		fprintf(candidatesfile, ",%f", *(candidate_dbl + i));
	}
	for(i = 0; i < numscenarios; i++){
		fprintf(candidatesfile, ",%f", *((*(candidate_UB_scen + candidate_ID)) + i));
	}
	fprintf(candidatesfile, ",%f\n", *(candidate_UB + candidate_ID));
	
	// Free space allocated to the scenario bounds
	free(*(candidate_UB_scen + candidate_ID));
	
	// Change candidate status
	*(candidate_status + candidate_ID) = CANDIDATE_EVALUATED;
	
}

/* Free all memory allocated to all candidates */
void free_memory_candidates(const int numcandidates, int *candidate_status,
	void **candidates, double **candidate_UB_scen)
{
	int i;
	for(i = 0; i < numcandidates; i++){
		free(*(candidates + i));
		if( *(candidate_status + i) == CANDIDATE_ACTIVE ){
			free(*(candidate_UB_scen + i));
		}
	}
}
