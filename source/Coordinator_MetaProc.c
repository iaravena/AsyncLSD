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

/* Allocate space in memory to hold registers for algorithm execution
 * and define pointers for easy access */
int allocate_memory_space_async_alg(const int numscenarios, const int numfirststagecols,
	const int initial_numcandidates, const int primal_queue_size,
	double **evaluation_timestamp,														// start time of last evaluation
	void **x_evaluated, double ***x_evaluated_scen_dbl,									// last evaluated multipliers
	void **v_evaluated, void ***v_evaluated_scen, double ***v_evaluated_scen_dbl,		// last evaluation results: LBs and subgradients
	void **u_evaluated, double **u_evaluated_dbl,										// last evaluated stoch program LB and f0_mu gradient (no time stamp needed)
	void **x_current, void ***x_current_scen, double ***x_current_scen_dbl,				// current multipliers
	void **x_nonscen, double **LBscen,
		double **msumx_evaluated, double **msumx_current, double **ucenter,				// last assembled non scenario task
	void ***candidates, int **candidate_status, double **candidate_UB,
		int **candidate_numevals, double ***candidate_UB_scen,							// candidate pointers (allocation is done dynamically as we go)
	int **primal_queue_candidate, int **primal_queue_scenario)							// primal queue, <candindate,scenario> pairs
{
	// Iterators
	int i;
	
	// Allocate buffers to store the data and communicate it
	*evaluation_timestamp = (double*) malloc( sizeof(double)*numscenarios );						// {timestampi : i=1,...,N}
	*x_evaluated = malloc( sizeof(double)*numfirststagecols*numscenarios );							// {xevali[] : i=1,...,N}
	*v_evaluated = malloc( (sizeof(int)*2 + sizeof(double)*(3+numfirststagecols))*numscenarios );	// {(sceni, XPRS status, worker time, fi_UB, fi_LB, vevali[]) : i=1,...,N}
	*u_evaluated = malloc( (sizeof(int)*2 + sizeof(double)*(3+numfirststagecols)) );				// (nonscen, XPRS status, worker time, f0_UB, f_LB, u[])
	*x_current = malloc( (sizeof(int) + sizeof(double)*numfirststagecols)*numscenarios );			// {(sceni, xcurrenti[]) : i=1,...,N}
	*x_nonscen = malloc( sizeof(int) + sizeof(double)*(1+3*numfirststagecols) );					// (NOSCENID, LBscen, msumx_evaluated[], msumx_current[], ucenter[])
	*candidate_status = (int*) malloc( sizeof(int)*initial_numcandidates );							// {statusj : j=1,...,numcandidates}
	*candidate_UB = (double*) malloc( sizeof(double)*initial_numcandidates );						// {UBj : j=1,...,numcandidates}
	*candidate_numevals = (int*) calloc(initial_numcandidates, sizeof(int));						// {num evaluated scenarios j : j=1,...,numcandidates}
	*primal_queue_candidate = (int*) malloc( sizeof(int)*primal_queue_size );						// {candidate}
	*primal_queue_scenario = (int*) malloc( sizeof(int)*primal_queue_size );						// {scenario}
	
	// Allocate pointers for easy access
	*x_evaluated_scen_dbl = (double**) malloc( sizeof(double*)*numscenarios );
	*v_evaluated_scen = (void**) malloc( sizeof(void*)*numscenarios );
	*v_evaluated_scen_dbl = (double**) malloc( sizeof(double*)*numscenarios );
	*x_current_scen = (void**) malloc( sizeof(void*)*numscenarios );
	*x_current_scen_dbl = (double**) malloc( sizeof(double*)*numscenarios );
	*candidates = (void**) malloc( sizeof(void*)*initial_numcandidates );
	*candidate_UB_scen = (double**) malloc( sizeof(double*)*initial_numcandidates );
	
	// Direct pointers to the corresponding arrays
	// x_eval
	**x_evaluated_scen_dbl = (double*) *x_evaluated;
	for(i = 1; i < numscenarios; i++){
		*((*x_evaluated_scen_dbl) + i) = (*((*x_evaluated_scen_dbl) + i - 1)) + numfirststagecols;
	}
	// v_eval, x_current
	**v_evaluated_scen = *v_evaluated;
	**x_current_scen = *x_current;
	for(i = 1; i < numscenarios; i++){
		*((*v_evaluated_scen) + i) = (void*) (((double*) (((int*) *((*v_evaluated_scen) + i - 1)) + 2)) + 3 + numfirststagecols);
		*((*x_current_scen) + i) = (void*) (((double*) (((int*) *((*x_current_scen) + i - 1)) + 1)) + numfirststagecols);
	}
	for(i = 0; i < numscenarios; i++){
		*((*v_evaluated_scen_dbl) + i) = ((double*) (((int*) *((*v_evaluated_scen) + i)) + 2)) + 3;
		*((*x_current_scen_dbl) + i) = (double*) (((int*) *((*x_current_scen) + i)) + 1);
	}
	// single element pointers
	*u_evaluated_dbl =  ((double*) (((int*) *u_evaluated) + 2)) + 3;
	*LBscen = (double*) (((int*) *x_nonscen) + 1);
	*msumx_evaluated = (*LBscen) + 1;
	*msumx_current = (*msumx_evaluated) + numfirststagecols;
	*ucenter = (*msumx_current) + numfirststagecols;
	
	// Initialize candidate status
	for(i = 0; i < initial_numcandidates; i++){
		*((*candidate_status) + i) = CANDIDATE_FREE_SLOT;
	}
	
	// Return success indicator
	return(0);
}

/* Free allocated memory for the asynchronous algorithm */
void free_memory_space_async_alg(double **evaluation_timestamp,
	void **x_evaluated, double ***x_evaluated_scen_dbl, void **v_evaluated, void ***v_evaluated_scen,
	double ***v_evaluated_scen_dbl, void **u_evaluated, void **x_current,
	void ***x_current_scen, double ***x_current_scen_dbl, void **x_nonscen,
	void ***candidates, int **candidate_status, double **candidate_UB,
		int **candidate_numevals, double ***candidate_UB_scen,
	int **primal_queue_candidate, int **primal_queue_scenario)
{
	// Issue free commands
	free(*evaluation_timestamp);
	free(*x_evaluated);
	free(*x_evaluated_scen_dbl);
	free(*v_evaluated);
	free(*v_evaluated_scen);
	free(*v_evaluated_scen_dbl);
	free(*u_evaluated);
	free(*x_current);
	free(*x_current_scen);
	free(*x_current_scen_dbl);
	free(*x_nonscen);
	free(*candidates);
	free(*candidate_status);
	free(*candidate_UB);
	free(*candidate_numevals);
	free(*candidate_UB_scen);
	free(*primal_queue_candidate);
	free(*primal_queue_scenario);
}

/* Initialization routine: compute sugradients and lower bound using period relaxation */
int initialize_subgradients_and_bounds(
	MPI_Datatype *dualscenjobmsg, MPI_Datatype *dualnonscenjobmsg, MPI_Datatype *dualjobresult,
	const int worldsize, struct slave *workers, MPI_Request *work_requests,
	const int k, const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, const double *update_scaling,
	void **x_current_scen, double **x_current_scen_dbl, double **x_evaluated_scen_dbl,
	void *x_nonscen, double *LBscen, double *msumx_evaluated, double *msumx_current, double *ucenter,
	double *evaluation_timestamp, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
	void *u_evaluated, double *u_evaluated_dbl, double *LBaux, double *LBk, double *gnorm2,
	int *numcandidates, int *maxnumcandidates, void ***candidates,
	int **candidate_status, double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
	int *num_arrivals, double *av_time_between_arrivals, double *stddev_time_between_arrivals,
	int *num_projections, double *av_time_projection, double *stddev_time_projection,
	struct timespec *startt, FILE *iterfile, const struct options *opt)
{
	// Iterators
	int i, j;
	
	// Create job queue
	struct generic_job *jobqueue = (struct generic_job*) malloc(sizeof(struct generic_job)*numscenarios);
	
	// Solve period relaxation in parallel using the workers
	struct timespec endt;
	MPI_Status result_status;
	void *result_buffer;
	double *arrival_time = (double*) malloc(sizeof(double)*numscenarios);
	
	// add dual jobs to the queue
	for(i = 0; i < numscenarios; i++){
		(*(jobqueue + i)).task = DUAL_SCEN_INIT;
		(*(jobqueue + i)).scenario = i;
	}
	
	// launch initial batch of jobs
	int worker_rank = 1;
	int queuepos = 0;
	while(queuepos < numscenarios && worker_rank < worldsize){
		// send job: scenario and multipliers
		if( MPI_Isend(*(x_current_scen + (*(jobqueue + queuepos)).scenario), 1, *dualscenjobmsg,
			worker_rank, (*(jobqueue + queuepos)).task, MPI_COMM_WORLD,
			work_requests + worker_rank) != MPI_SUCCESS ){
			printf("\nError while sending initialization job through MPI_Isend.\n");
			return(-1);
		}
		// keep track of what each worker is doing
		current_utc_time(&endt);
		(*(workers + worker_rank)).status = WORKER_BUSY_DUAL;
		(*(workers + worker_rank)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
					((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
		(*(workers + worker_rank)).current_job.task = (*(jobqueue + queuepos)).task;
		(*(workers + worker_rank)).current_job.scenario = (*(jobqueue + queuepos)).scenario;
		// increase count
		worker_rank++;
		queuepos++;
	}
	
	// collect results and execute remaining jobs (if any)
	int collected_results = 0;
	while( collected_results < numscenarios ){
		// get next result type
		if( MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&result_status) != MPI_SUCCESS ){
			printf("\nError while probing for result message (MPI_Probe).\n");
			return(-1);
		}
		// check that result type corresponds to initializations
		if( result_status.MPI_TAG != RESULT_DUAL_SCEN_INIT ){
			printf("\nUnexpected tag %d received.\n", result_status.MPI_TAG);
			return(-1);
		}
		// receive result
		result_buffer = *(v_evaluated_scen + (*(workers + result_status.MPI_SOURCE)).current_job.scenario);
		if( MPI_Recv(result_buffer, 1, *dualjobresult,
			result_status.MPI_SOURCE, result_status.MPI_TAG,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
			printf("\nError while receiving message %d from worker %d.\n",
				result_status.MPI_TAG, result_status.MPI_SOURCE);
			return(-1);
		}
		// print received message
		current_utc_time(&endt);
		if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
			printf("%5.1f\t%d\t%5.1f\t%d\tINIT\t%d\t-----\t%5.2f\t-----\t-----\t-----\n",
				((double) (endt.tv_sec - (*startt).tv_sec)) +
					((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
				result_status.MPI_SOURCE,
				*((double*) (((int*) result_buffer) + 2)),
				k,
				*((int*) result_buffer),
				*(((double*) (((int*) result_buffer) + 2)) + 2)/1E6);
		// write iteration data
		fprintf(iterfile, "%f,%d,%f,%d,INIT,%d,%f,%f,,\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			result_status.MPI_SOURCE,
			*((double*) (((int*) result_buffer) + 2)),
			k,
			*((int*) result_buffer),
			*(((double*) (((int*) result_buffer) + 2)) + 1)/1E6,
			*(((double*) (((int*) result_buffer) + 2)) + 2)/1E6);
		// copy evaluated x into the corresponding register
		*(evaluation_timestamp + (*(workers + result_status.MPI_SOURCE)).current_job.scenario) =
			(*(workers + result_status.MPI_SOURCE)).dispatch_timestamp;
		memcpy((void*) *(x_evaluated_scen_dbl + (*(workers + result_status.MPI_SOURCE)).current_job.scenario),
			(void*) *(x_current_scen_dbl + (*(workers + result_status.MPI_SOURCE)).current_job.scenario),
			sizeof(double)*numfirststagecols);
		// register arrival time
		*(arrival_time + collected_results) = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
		// mark worker as free
		(*(workers + result_status.MPI_SOURCE)).status = WORKER_FREE;
		// start next process in the queue if necessary
		if(queuepos < numscenarios){
			// send job: scenario and multipliers
			if( MPI_Isend(*(x_current_scen + (*(jobqueue + queuepos)).scenario), 1, *dualscenjobmsg,
				result_status.MPI_SOURCE, (*(jobqueue + queuepos)).task, MPI_COMM_WORLD,
				work_requests + result_status.MPI_SOURCE) != MPI_SUCCESS ){
				printf("\nError while sending initialization job through MPI_Isend.\n");
				return(-1);
			}
			// keep track of what each worker is doing
			current_utc_time(&endt);
			(*(workers + result_status.MPI_SOURCE)).status = WORKER_BUSY_DUAL;
			(*(workers + result_status.MPI_SOURCE)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
					((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
			(*(workers + result_status.MPI_SOURCE)).current_job.task = (*(jobqueue + queuepos)).task;
			(*(workers + result_status.MPI_SOURCE)).current_job.scenario = (*(jobqueue + queuepos)).scenario;
			// increase count
			queuepos++;
		}
		// increment count
		collected_results++;
	}
	// free queue and flush registers
	free(jobqueue);
	fflush(stdout);
	fflush(iterfile);
	
	// compute arrival time statitics and free allocated vector
	double sumt = 0, sumt2 = 0;
	for(i = 1; i < numscenarios; i++){
		sumt += *(arrival_time + i) - *(arrival_time + i - 1);
		sumt2 += pow(*(arrival_time + i) - *(arrival_time + i - 1), 2.0);
	}
	*num_arrivals = numscenarios - 1;
	*av_time_between_arrivals = sumt/(*num_arrivals);
	*stddev_time_between_arrivals = pow((sumt2 - (*num_arrivals)*pow(*av_time_between_arrivals, 2.0))/((double) (*num_arrivals)-1.0), 0.5);
	free(arrival_time);
	
	// Launch non-scenario problem to compute a lower bound on the stochastic program
	
	// update x_nonscen
	*LBscen = 0.0;
	for(i = 0; i < numscenarios; i++)
		*LBscen += *(scenprobabilities + i) * *(((double*) (((int*) *(v_evaluated_scen + i)) + 2)) + 2);
	for(j = 0; j < numfirststagecols; j++){
		*(msumx_evaluated + j) = 0.0;
		*(ucenter + j) = 0.0;
	}
	for(i = 0; i < numscenarios; i++){
		for(j = 0; j < numfirststagecols; j++){
			*(msumx_evaluated + j) -= *(scenprobabilities + i) * *((*(x_evaluated_scen_dbl + i)) + j);
			*(ucenter + j) += *(scenprobabilities + i) * *((*(v_evaluated_scen_dbl + i)) + j);
		}
	}
	
	// launch non-scenario job
	const int worker_rank_nonscen = 1;							// i.e. launching using first worker
	if( MPI_Isend(x_nonscen, 1, *dualnonscenjobmsg,
		worker_rank_nonscen, DUAL_NONSCEN_MILP, MPI_COMM_WORLD,
		work_requests + worker_rank_nonscen) != MPI_SUCCESS ){
		printf("\nError while sending dual job through MPI_Isend.\n");
		return(-1);
	}
	
	// keep track of what each worker is doing
	current_utc_time(&endt);
	(*(workers + worker_rank_nonscen)).status = WORKER_BUSY_DUAL;
	(*(workers + worker_rank_nonscen)).dispatch_timestamp = ((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC;
	(*(workers + worker_rank_nonscen)).current_job.task = DUAL_NONSCEN_MILP;
	(*(workers + worker_rank_nonscen)).current_job.scenario = NOSCENID;
	
	// Create initial set of primal candidates
	//double average_init_primals = (worldsize-1)*(1-(*opt).dual_share_start);						// binomial mean
	//double sd_init_primals = (worldsize-1)*(1-(*opt).dual_share_start)*(*opt).dual_share_start;	// binomial variance
	//double primal_queue_critical_lenght = (int) ceil((1-(*opt).dual_share_start)*(worldsize-1));	// critical lenght, if there are less primal jobs queued we will add new jobs
	//int numinitprimal = (int) ceil((average_init_primals + sd_init_primals +
	//	primal_queue_critical_lenght)/((double) numscenarios));		// mu + sigma + critical length: this gives a 84% chance that I will have enough jobs in the first round
	if( (*opt).primal_recovery == PRIMAL_RECOVERY_FIFO
		|| (*opt).primal_recovery == PRIMAL_RECOVERY_RND
		|| (*opt).primal_recovery == PRIMAL_RECOVERY_LIFO ){
		// add all possible candidates from the initialization results
		int numinitprimal = (int) ceil(2*ceil((worldsize-1)*(1.0-(*opt).dual_share_start))/((double) numscenarios));
		int scen;
		for(scen = 0; scen < numscenarios; scen++){
			if( add_primal_candidate(numscenarios, numfirststagecols, *(v_evaluated_scen_dbl + scen),
				numcandidates, maxnumcandidates, candidates, candidate_status,
				candidate_UB, candidate_numevals, candidate_UB_scen, opt) < 0 ){
				printf("\nError adding a new primal candidate into the list.\n");
				return(-1);
			}
		}
		// check number of added candidates
		if( *numcandidates < numinitprimal + 1 ){
			printf("\nFailed to generate %d candidates during initializations."
			"\nMaybe too many processors were launched for this problem.\n", numinitprimal + 1);
			return(-1);
		}
		/*
		// activate initial primal candidates
		*primal_queue_numjobs = 0;
		for(i = 0; i < numinitprimal; i++){
			if( activate_primal_candidate(i, numscenarios, *numcandidates,
				*candidate_status, *candidate_UB_scen, primal_queue_size,
				primal_queue_numjobs, primal_queue_first, primal_queue_last,
				primal_queue_candidate, primal_queue_scenario) != 0 ){
				printf("\nError while trying to activate candidate %d.\n", i);
				return(-1);
			}
		}
		*/
	}else if( (*opt).primal_recovery == PRIMAL_RECOVERY_IS ){
		;
		/*
		// get numinitprimal + 2 candidates using the IS heuristic
		int n = 0;
		int worker_rank_project = worldsize-1;
		void *new_candidate = malloc(sizeof(int)*2 + sizeof(double)*numfirststagecols);		// (candj, scenario dummy, uj[])
		double *new_candidate_dbl = (double*) (((int*) new_candidate) + 2);
		void *primal_proj_receiver = malloc( sizeof(int) + sizeof(double)*(1+numfirststagecols) );
		double *new_candidate_feasible = ((double*) (((int*) primal_proj_receiver) + 1)) + 1;
		double sumtproj = 0, sumt2proj = 0;
		while( *numcandidates < numinitprimal + 1 && n < MAX_NTRIALS_PRIMAL_GEN ){
			// get infeasible new primal
			generate_new_IS_primal_candidate(numscenarios, scenprobabilities, numfirststagecols,
				v_evaluated_scen, v_evaluated_scen_dbl, new_candidate_dbl, opt);
			// send and receive projection task
			*((int*) new_candidate) = NOCANDID;
			*(((int*) new_candidate) + 1) = NOSCENID;
			if( MPI_Send(new_candidate, 1, *primaljobmsg,
				worker_rank_project, PRIMAL_FEAS_PROJ, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending projection job through MPI_Send.\n");
				return(-1);
			}
			if( MPI_Recv(primal_proj_receiver, 1, *primalprojjobresult,
				worker_rank_project, RESULT_PRIMAL_FEAS_PROJ,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from worker %d.\n",
					RESULT_PRIMAL_FEAS_PROJ, worker_rank_project);
				return(-1);
			}
			// print received message
			current_utc_time(&endt);
			if((*opt).verbosity >= VERB_COORDINATOR_FUNCS)
				printf("%5.1f\t%d\t%5.1f\t%d\tPROJ\t-----\t-----\t-----\t-----\t-----\t-----\n",
					((double) (endt.tv_sec - (*startt).tv_sec)) +
						((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
					worker_rank_project, *((double*) (((int*) primal_proj_receiver) + 1)), k);
			// write iteration data
			fprintf(iterfile, "%f,%d,%f,%d,PROJ,,,,,\n",
				((double) (endt.tv_sec - (*startt).tv_sec)) +
					((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
				worker_rank_project, *((double*) (((int*) primal_proj_receiver) + 1)), k);
			// add primal to the list
			if( add_primal_candidate(numscenarios, numfirststagecols, new_candidate_feasible,
				numcandidates, maxnumcandidates, candidates, candidate_status,
				candidate_UB, candidate_numevals, candidate_UB_scen, opt) < 0 ){
				printf("\nError adding a new primal candidate into the list.\n");
				return(-1);
			}
			// increase attempts counter and register solution time
			n++;
			sumtproj += *((double*) (((int*) primal_proj_receiver) + 1));
			sumt2proj += pow(*((double*) (((int*) primal_proj_receiver) + 1)), 2.0);
		}
		// compute projection time statistics
		*num_projections = n;
		*av_time_projection = sumtproj/((double) n);
		*stddev_time_projection = pow((sumt2proj - n*pow(*av_time_projection, 2.0))/((double) n - 1.0), 0.5);
		// free space allocated to the candidate
		free(new_candidate);
		free(primal_proj_receiver);
		*/
	}else{
		printf("\nUnrecognized primal recovery method %d.\n", (*opt).primal_recovery);
		return(-1);
	}
	
	// Receive results of non-scenario problem, update lower bound and subgradient norm
	if( MPI_Recv(u_evaluated, 1, *dualjobresult,
		worker_rank_nonscen, RESULT_DUAL_NONSCEN_MILP,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
		printf("\nError while receiving message %d from worker %d.\n",
			RESULT_DUAL_NONSCEN_MILP, worker_rank_nonscen);
		return(-1);
	}
	
	// release worker
	(*(workers + worker_rank_nonscen)).status = WORKER_FREE;
	
	// compute bounds and subgradient norm
	struct timespec intstartt;
	*LBaux = *(((double*) (((int*) u_evaluated) + 2)) + 2);
	*LBk = fmax(*LBk, *LBaux);
	current_utc_time(&intstartt);
	*gnorm2 = subgradnorm2(numscenarios, scenprobabilities, numfirststagecols, update_scaling,
		v_evaluated_scen_dbl, u_evaluated_dbl);
	current_utc_time(&endt);
	double gnorm2time = 1000*(((double) (endt.tv_sec - intstartt.tv_sec)) +
				((double) (endt.tv_nsec - intstartt.tv_nsec))*NANOSEC2SEC);
	
	// print current state of problem
	current_utc_time(&endt);
	if((*opt).verbosity >= VERB_COORDINATOR_ONLY)
		printf("%5.1f\t%d\t%5.1f\t%d\tINIT\tAll\t-----\t%5.2f\t-----%%\t%5.1f\t%5.3fm\n",
			((double) (endt.tv_sec - (*startt).tv_sec)) +
				((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
			worker_rank_nonscen, *((double*) (((int*) u_evaluated) + 2)),
			k, *LBaux/1E6, *gnorm2, gnorm2time);
	fflush(stdout);
	// write iteration data
	fprintf(iterfile, "%f,%d,%f,%d,INIT,All,,%f,%f,\n",
		((double) (endt.tv_sec - (*startt).tv_sec)) +
			((double) (endt.tv_nsec - (*startt).tv_nsec))*NANOSEC2SEC,
		worker_rank_nonscen, *((double*) (((int*) u_evaluated) + 2)),
		k, *LBaux/1E6, *LBk/1E6);
	fflush(iterfile);
	
	// Return success indicator
	return(0);
}
