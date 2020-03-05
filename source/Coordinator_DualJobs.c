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

/* Update multipliers and launch dual scenario subproblem */
int update_and_launch_dual_scen_subprob(const int k, const int scenario,
	const double gnorm2, const double LB, const double UB,
	const int numscenarios, const double *scenprobabilities,							// required for stepsize computation
	const int numfirststagecols, const double *update_scaling,
	void **x_current_scen, double **x_current_scen_dbl,
	double *msumx_current, double **v_evaluated_scen_dbl, double *u_evaluated_dbl,		// required for updating
	const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *dualscenjobmsg, struct timespec *startt, const struct options *opt)
{
	// Check whether the assigned worker is actually free
	if( (*(workers + assigned_worker)).status != WORKER_FREE ){
		printf("\nWorker assigned for next dual iteration is not free."
			"\nWorker status: %d (it should be %d)\n",
			(*(workers + assigned_worker)).status, WORKER_FREE);
		return(-1);
	}
	
	// Select scenario at random and perform update
	int update_scen;
	if(scenario == RANDOM_SCENARIO){
		update_scen = rand() % numscenarios;
	}else{
		update_scen = scenario;
	}
	double lambdapi = stepsize(k, gnorm2, LB, UB, numscenarios, opt) * *(scenprobabilities + update_scen);
	double delta_x;
	int j;
	for(j = 0; j < numfirststagecols; j++){
		if( *(update_scaling + j) > NO_UPDATE_SCALING ){
			delta_x = lambdapi * (*(update_scaling + j)) *
				( (*((*(v_evaluated_scen_dbl + update_scen)) + j)) - *(u_evaluated_dbl + j) );
		}else{
			delta_x = lambdapi * ( (*((*(v_evaluated_scen_dbl + update_scen)) + j)) - *(u_evaluated_dbl + j) );
		}
		*((*(x_current_scen_dbl + update_scen)) + j) += delta_x;
		*(msumx_current + j) -= *(scenprobabilities + update_scen) * delta_x;
	}
	
	// Allocate task details, including multipliers, in the workers object
	struct timespec endt;
	current_utc_time(&endt);
	(*(workers + assigned_worker)).status = WORKER_BUSY_DUAL;
	(*(workers + assigned_worker)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
	(*(workers + assigned_worker)).current_job.task = DUAL_SCEN_MILP;
	(*(workers + assigned_worker)).current_job.scenario = update_scen;
	(*(workers + assigned_worker)).current_job.x_task = malloc( sizeof(int) + sizeof(double)*numfirststagecols );
	memcpy((*(workers + assigned_worker)).current_job.x_task,
		(void*) *(x_current_scen + update_scen), sizeof(int) + sizeof(double)*numfirststagecols);
	
	// Launch evaluation of updated scenario (note that the task is sent using non-blocking communications)
	if( MPI_Isend((*(workers + assigned_worker)).current_job.x_task, 1, *dualscenjobmsg,
		assigned_worker, (*(workers + assigned_worker)).current_job.task, MPI_COMM_WORLD,
		work_requests + assigned_worker) != MPI_SUCCESS ){
		printf("\nError while sending dual scenario subproblem job through MPI_Isend.\n");
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Posprocess dual scenario subproblem */
int posprocess_dual_scen_subprob(const int k, const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, double *evaluation_timestamp, double **x_evaluated_scen_dbl,
	double *LBscen, double *msumx_evaluated, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
	double *ucenter, MPI_Status *result_status, struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *dualjobresult, int *numcandidates, int *maxnumcandidates, void ***candidates,
	int **candidate_status, double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
	struct timespec *startt, FILE *iterfile, const struct options *opt)
{
	// Check whether we will receive updated information and define receiver accordingly
	void *result_buffer;
	int incumbent_worker = (*result_status).MPI_SOURCE;
	int incumbent_scenario = (*(workers + incumbent_worker)).current_job.scenario;
	double old_LB_evaluated;
	int newinfo = (*(workers + incumbent_worker)).dispatch_timestamp > *(evaluation_timestamp + incumbent_scenario);
	if( newinfo ){
		// information is new! -> update v_evaluated and x_evaluated registers
		int i;
		double *old_x_evaluated = *(x_evaluated_scen_dbl + incumbent_scenario);
		double *new_x_evaluated = (double*) (((int*) (*(workers + incumbent_worker)).current_job.x_task) + 1);
		for(i = 0; i < numfirststagecols; i++){
			*(msumx_evaluated + i) -= *(scenprobabilities + incumbent_scenario) *
				((*(new_x_evaluated + i)) - (*(old_x_evaluated + i)));
		}
		memcpy((void*) old_x_evaluated, (void*) new_x_evaluated, sizeof(double)*numfirststagecols);
		*(evaluation_timestamp + incumbent_scenario) = 
			(*(workers + incumbent_worker)).dispatch_timestamp;
		// substract scenario component from center
		if( (*opt).nonscen_center_type == CENTER_AVERAGE_SUBGRADIENT ){
			for(i = 0; i < numfirststagecols; i++){
				*(ucenter + i) -= *(scenprobabilities + incumbent_scenario) *
					*((*(v_evaluated_scen_dbl + incumbent_scenario)) + i);
			}
		}
		// prepare buffer for reception
		result_buffer = *(v_evaluated_scen + incumbent_scenario);
		old_LB_evaluated = *(((double*) (((int*) result_buffer) + 2)) + 2);
	}else{
		// information is older than what we already have
		// -> receive in internal buffer
		result_buffer = malloc(sizeof(int)*2 + sizeof(double)*(3+numfirststagecols));
	}
	
	// Receive the message
	if( MPI_Recv(result_buffer, 1, *dualjobresult,
		(*result_status).MPI_SOURCE, (*result_status).MPI_TAG,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
		printf("\nError while receiving message %d from worker %d.\n",
			(*result_status).MPI_TAG, (*result_status).MPI_SOURCE);
		return(-1);
	}
	
	// Check that the received scenario coincides with the one registered by the worker
	if(incumbent_scenario != *((int*) result_buffer)){
		printf("\nInconsistent reception of scenarios. Expected scenario %d, received scenario %d.\n",
			incumbent_scenario, *((int*) result_buffer));
		return(-1);
	}
	
	// Print received message
	struct timespec endt;
	current_utc_time(&endt);
	if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
		printf("%5.1f\t%d\t%5.1f\t%d\tDUAL\t%d\t%5.2f\t%5.2f\t-----\t-----\t-----\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			incumbent_worker,
			*((double*) (((int*) result_buffer) + 2)),
			k,
			*((int*) result_buffer),
			*(((double*) (((int*) result_buffer) + 2)) + 1)/1E6,
			*(((double*) (((int*) result_buffer) + 2)) + 2)/1E6);
	// write iteration data
	fprintf(iterfile, "%f,%d,%f,%d,DUAL,%d,%f,%f,,\n",
		((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
		incumbent_worker,
		*((double*) (((int*) result_buffer) + 2)),
		k,
		*((int*) result_buffer),
		*(((double*) (((int*) result_buffer) + 2)) + 1)/1E6,
		*(((double*) (((int*) result_buffer) + 2)) + 2)/1E6);
	
	// Complete update of LBscen
	if( newinfo ){
		*LBscen += *(scenprobabilities + incumbent_scenario) *
			((*(((double*) (((int*) result_buffer) + 2)) + 2)) - old_LB_evaluated);
	}
	
	// Complete update of ucenter
	if( newinfo && (*opt).nonscen_center_type == CENTER_AVERAGE_SUBGRADIENT ){
		int i;
		for(i = 0; i < numfirststagecols; i++){
			*(ucenter + i) += *(scenprobabilities + incumbent_scenario) *
				*((*(v_evaluated_scen_dbl + incumbent_scenario)) + i);
		}
	}
	
	// Add candidate to the list
	if( (*opt).primal_recovery == PRIMAL_RECOVERY_FIFO
		|| (*opt).primal_recovery == PRIMAL_RECOVERY_RND
		|| (*opt).primal_recovery == PRIMAL_RECOVERY_LIFO ){
		//printf("\nCandidated. Before: %d", *numcandidates);
		if( add_primal_candidate(numscenarios, numfirststagecols,
			((double*) (((int*) result_buffer) + 2)) + 3,
			numcandidates, maxnumcandidates, candidates, candidate_status,
			candidate_UB, candidate_numevals, candidate_UB_scen, opt) < 0){
			printf("\nError adding a new primal candidate into the list.\n");
			return(-1);
		}
		//printf("\tAfter: %d\n", *numcandidates);
	}
	
	// Release worker
	(*(workers + incumbent_worker)).status = WORKER_FREE;
	free((*(workers + incumbent_worker)).current_job.x_task);	// here is where multipliers were stored
	
	// Release internal memory
	if( !newinfo ){
		free(result_buffer);
	}
	
	// Return success indicator
	return(0);
}

/* Launch non-scenario dual subproblems */
int launch_dual_nonscen_subprob(const int numfirststagecols, const void *x_nonscen,
	const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
	MPI_Datatype *dualnonscenjobmsg, struct timespec *startt, const struct options *opt)
{
	// Check whether the assigned worker is actually free
	if( (*(workers + assigned_worker)).status != WORKER_FREE ){
		printf("\nWorker assigned for solving dual non-scenario subproblems is not free."
			"\nWorker status: %d (it should be %d)\n",
			(*(workers + assigned_worker)).status, WORKER_FREE);
		return(-1);
	}
	
	// Allocate task details, including multipliers in the workers object
	struct timespec endt;
	current_utc_time(&endt);
	(*(workers + assigned_worker)).status = WORKER_BUSY_DUAL;
	(*(workers + assigned_worker)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
	(*(workers + assigned_worker)).current_job.task = DUAL_NONSCEN_MILP;
	(*(workers + assigned_worker)).current_job.scenario = NOSCENID;
	(*(workers + assigned_worker)).current_job.x_task = malloc( sizeof(int) + sizeof(double)*(1+3*numfirststagecols) );
	memcpy((*(workers + assigned_worker)).current_job.x_task, x_nonscen,
		sizeof(int) + sizeof(double)*(1+3*numfirststagecols) );
	
	// launch non-scenario job
	if( MPI_Isend((*(workers + assigned_worker)).current_job.x_task, 1,
		*dualnonscenjobmsg, assigned_worker,
		(*(workers + assigned_worker)).current_job.task, MPI_COMM_WORLD,
		work_requests + assigned_worker) != MPI_SUCCESS ){
		printf("\nError while sending dual non-scenario job through MPI_Isend.\n");
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Posprocess dual non-scenario subproblem */
int posprocess_dual_nonscen_subprob(const int k, const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, const double *update_scaling,
	void *u_evaluated, double *u_evaluated_dbl, double **v_evaluated_scen_dbl,
	const double UBk, double *LBaux, double *LBk, double *gnorm2, MPI_Status *result_status,
	struct slave *workers, MPI_Request *work_requests, MPI_Datatype *dualjobresult,
	struct timespec *startt, FILE *iterfile, const struct options *opt)
{
	// Receive results of non-scenario problem
	int incumbent_worker = (*result_status).MPI_SOURCE;
	if( MPI_Recv(u_evaluated, 1, *dualjobresult,
		incumbent_worker, RESULT_DUAL_NONSCEN_MILP,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
		printf("\nError while receiving message %d from worker %d.\n",
			RESULT_DUAL_NONSCEN_MILP, incumbent_worker);
		return(-1);
	}
	
	// Check that the received scenario coincides with NOSCENID
	if(NOSCENID != *((int*) u_evaluated)){
		printf("\nInconsistent reception of subproblem solution. Expected scenario %d, received scenario %d.\n",
			NOSCENID, *((int*) u_evaluated));
		return(-1);
	}
	
	// compute bounds and subgradient norm
	struct timespec intstartt, endt;
	*LBaux = *(((double*) (((int*) u_evaluated) + 2)) + 2);
	*LBk = fmax(*LBk, *LBaux);
	current_utc_time(&intstartt);
	*gnorm2 = (1-(*opt).dual_gnorm2_filter)*subgradnorm2(numscenarios,
		scenprobabilities, numfirststagecols, update_scaling, v_evaluated_scen_dbl, u_evaluated_dbl) + 
		(*opt).dual_gnorm2_filter*(*gnorm2);
	current_utc_time(&endt);
		double gnorm2time = 1000*(((double) (endt.tv_sec - intstartt.tv_sec)) +
				((double) (endt.tv_nsec - intstartt.tv_nsec))*NANOSEC2SEC);
	
	// print current state of problem
	double dgap = 100*(UBk - *LBk)/(1E-10 + fmax(fabs(UBk),fabs(*LBk)));
	current_utc_time(&endt);
	if((*opt).verbosity >= VERB_COORDINATOR_ONLY)
		printf("%5.1f\t%d\t%5.1f\t%d\tDUAL\tAll\t%5.2f\t%5.2f\t%5.2f%%\t%5.1f\t%5.3fm\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			incumbent_worker, *((double*) (((int*) u_evaluated) + 2)),
			k, UBk/1E6, *LBaux/1E6, dgap, *gnorm2, gnorm2time);
	fflush(stdout);
	// write iteration data
	fprintf(iterfile, "%f,%d,%f,%d,DUAL,All,%f,%f,%f,%f\n",
		((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
		incumbent_worker, *((double*) (((int*) u_evaluated) + 2)),
		k, UBk/1E6, *LBaux/1E6, *LBk/1E6, dgap);
	fflush(iterfile);
	
	// release worker
	(*(workers + incumbent_worker)).status = WORKER_FREE;
	free((*(workers + incumbent_worker)).current_job.x_task);
	
	// Return success indicator
	return(0);
}