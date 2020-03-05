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
 * Worker routines derivated from main
 */

/* Asynchronous algorithm header */
#include "AsyncHeader.h"

/* Routines to solve subproblems */
int solve_init_scenario_subproblem(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, const double *dual_multipliers,
	int *probstatus, double *UB, double *LB, double *firststagesolution, const struct options* global_opt);
int solve_dual_scenario_subproblem(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int numfirststagecols, const int *firststagecols, const double *firststagelb, const double *firststageub,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, const double *dual_multipliers,
	int *probstatus, double *UB, double *LB, double *firststagesolution, const struct options* global_opt);
int solve_nonscenario_subproblems(const char *workdir, const struct string_buffer *colnames,
	const struct string_buffer *rownames, const int *rowtstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	XPRSprob *workingproblem, const double *msum_delayed_multipliers,
	const double *msum_current_multipliers, const double *ucenter,
	int *probstatus, double *UB, double *LB, double *solution, const struct options *global_opt);
int project_infeasible_candidate(const char *workdir, const struct string_buffer *colnames,
	const struct string_buffer *rownames, const int *rowtstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	XPRSprob *workingproblem, const double *infeasible_candidate,
	int *probstatus, double *solution, const struct options* global_opt);
int evaluate_primal_candidate(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int numfirststagecols, const int *firststagecols,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, double *first_stage_candidate,
	int *probstatus, double *UB, double *LB, const struct options* global_opt);
void export_troublemaker(const char *workdir, const char *filename, const char *format, XPRSprob *problem);

/* Scenario worker */
int ScenarioWorker(const int worldrank, int *argc, char ***argv)
{
	// Time registers
	struct timespec startt, endt;
	
	// Options for execution
	struct options global_opt;
	if( read_options(*((*argv) + 4), &global_opt) != 0 ){
		printf("\nError while reading options. Program will exit.\n");
		return(-1);
	}
	
	// Synchronization point: waiting for files to be ready
	MPI_Barrier( MPI_COMM_WORLD );
	if(global_opt.verbosity >= VERB_WORKER)
		printf("\nP%d: Files ready for reading.", worldrank);
	
	/*
	// debug only
	MPI_Barrier( MPI_COMM_WORLD );
	// debug only
	*/
	
	// Initialize Xpress-Optimizer
	char *message = (char*) malloc(sizeof(char)*512);
	if( XPRSinit(NULL) ) {
		XPRSgetlicerrmsg(message,512);
		printf("%s\n", message);
		return(-1);
	}
	free(message);
	
	// Variables and pointers for manipulating subproblems
	struct string_buffer row_names, col_names, time_stages, scenarios, scens_parent;
	int *rows_tstage, *cols_tstage, *scens_branch_tstage, *scens_start_line, *scens_end_line,
		firststageind, numfirststagecols, *firststagecols;
	double *scens_probability, *firststagelb, *firststageub;
	
	// Main problem of the worker
	XPRSprob prob;
	XPRScreateprob(&prob);
	XPRSsetintcontrol(prob, XPRS_SCALING, 0);			// avoids automatic scaling when reading the problem
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		XPRSaddcbmessage(prob, Message, NULL, 0);		// defines callback function for printing log from Xpress (usefull only in Windows)
	#endif
	
	// Load Xpress parameters
	{
		int i;
		for(i = 0; i < global_opt.scenprob_options.num_intoptions; i++){
			//printf("\nSetting %d to %d",
			//	*(global_opt.scenprob_options.intoptions + i),
			//	*(global_opt.scenprob_options.intoptions_value + i));
			//fflush(stdout);
			XPRSsetintcontrol(prob, *(global_opt.scenprob_options.intoptions + i),
				*(global_opt.scenprob_options.intoptions_value + i));
		}
		for(i = 0; i < global_opt.scenprob_options.num_dbloptions; i++){
			XPRSsetdblcontrol(prob, *(global_opt.scenprob_options.dbloptions + i),
				*(global_opt.scenprob_options.dbloptions_value + i));
		}
	}
	
	// Override verbosity level
	if(global_opt.verbosity > VERB_WORKER){
		XPRSsetintcontrol(prob, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(prob, XPRS_OUTPUTLOG, 4);
	}
	
	// Read problem from SMPS files
	if(global_opt.verbosity >= VERB_WORKER){
		printf("\nP%d: Reading SMPS instance into memory...\n", worldrank);
		current_utc_time(&startt);
	}
	if( read_SMPS(*((*argv) + 1), *((*argv) + 3),
		&prob, &row_names, &col_names,
		&time_stages, &rows_tstage, &cols_tstage,
		&scenarios, &scens_parent, &scens_branch_tstage,
		&scens_start_line, &scens_end_line, &scens_probability,
		&firststageind, &numfirststagecols, &firststagecols,
		&firststagelb, &firststageub, &global_opt) != 0 ) {
		printf("\nError while reading SMPS files.\n");
		return(-1);
	}
	if(global_opt.verbosity >= VERB_WORKER){
		int nonzeros, intvars;
		XPRSgetintattrib(prob, XPRS_ELEMS, &nonzeros);
		XPRSgetintattrib(prob, XPRS_MIPENTS, &intvars);
		current_utc_time(&endt);
		printf("\nP%d: CORE, STOCH and TIME files read in %.2f secs."
			"\nP%d: CORE problem has %d rows, %d cols (%d non-zeros) and %d integer variables."
			"\nP%d: %d scenarios and %d time stages detected in STOCH and TIME files.\n",
			worldrank, ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC,
			worldrank, row_names.num_elems, col_names.num_elems, nonzeros, intvars,
			worldrank, scenarios.num_elems, time_stages.num_elems);
	}
	
	// Synchronization point: waiting for files to be read
	MPI_Barrier( MPI_COMM_WORLD );
	if(global_opt.verbosity >= VERB_WORKER)
		printf("\nP%d: All files read. Proceeding with algorithm execution...", worldrank);
	
	// Define messages
	MPI_Datatype dualscenjobmsg, dualnonscenjobmsg, dualjobresult,
		primaljobmsg, primalprojjobresult, primaljobresult;
	if( create_message_types(numfirststagecols,
		&dualscenjobmsg, &dualnonscenjobmsg, &dualjobresult,
		&primaljobmsg, &primalprojjobresult, &primaljobresult) != 0 ){
		printf("\nError while creating MPI message types.\n");
		return(-1);
	}
	
	// Waiting for a task
	MPI_Status task_status;
	void *task_dual_scen_buffer = malloc( sizeof(int) + sizeof(double)*numfirststagecols );				// (scen, xi[])
	void *task_dual_nonscen_buffer = malloc( sizeof(int) + sizeof(double)*(1 + 3*numfirststagecols) );	// (nonscen, LBscen^delayed, -sum_i x^delayed_i[], -sum_i x^current_i[], incumbent)
	void *task_primal_buffer = malloc( sizeof(int)*2 + sizeof(double)*numfirststagecols );				// (candID, scen, candidate[])
	void *dual_result_buffer = malloc( sizeof(int)*2 + sizeof(double)*(3+numfirststagecols) );			// (scen, XPRS status, worker time, fi_UB, fi_LB, vi[]) or (nonscen, XPRS status, f0_UB, f0_LB, u[])
	void *primal_proj_result_buffer = malloc( sizeof(int) + sizeof(double)*(1+numfirststagecols) );		// (candID, worker time, feasible candidate[])
	void *primal_result_buffer = malloc( sizeof(int)*3 + sizeof(double)*3 );							// (candID, sceni, XPRS status, worker time, fi(cand)_UB, fi(cand)_LB)
	int jobcounter = 0; 
	while( 1 ){
		
		// probe next task
		if( MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &task_status) != MPI_SUCCESS ){
			printf("\nError while probing for task message (MPI_Probe).\n");
			return(-1);
		}
		
		jobcounter++;
		if(global_opt.verbosity >= VERB_WORKER){
			printf("\nP%d: Received its task No. %d. TAG: %d\n", worldrank, jobcounter, task_status.MPI_TAG);
			fflush(stdout);
			}
		
		// start counting worker time
		current_utc_time(&startt);
		
		// what to do next depends on the TAG ...
		if( task_status.MPI_TAG == DUAL_SCEN_INIT ){		// dual initialization scenario subproblem
			
			// receive scenario and multipliers
			if( MPI_Recv(task_dual_scen_buffer, 1, dualscenjobmsg, 0, task_status.MPI_TAG,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from Coordinator.\n", task_status.MPI_TAG);
				return(-1);
			}
			
			// put scenario into result buffer
			*((int*) dual_result_buffer) = *((int*) task_dual_scen_buffer);
			
			// solve period relaxation with the provided multipliers
			if(global_opt.verbosity >= VERB_WORKER)
				printf("\nP%d: Solving period relaxation of dual scenario subproblem for scenario %s.\n",
					worldrank, scenarios.buffer + (scenarios.atom_len+1) * *((int*) task_dual_scen_buffer));
			if( solve_init_scenario_subproblem(*(*argv + 1), *(*argv + 3),
				&row_names, &col_names,
				&time_stages, rows_tstage, cols_tstage,
				firststageind, numfirststagecols, firststagecols,
				firststagelb, firststageub,
				&scenarios, &scens_parent, scens_branch_tstage,
				scens_start_line, scens_end_line, scens_probability,
				&prob,
				*((int*) task_dual_scen_buffer),												// scenario number
				(double*) (((int*) task_dual_scen_buffer) + 1),									// multipliers
				((int*) dual_result_buffer) + 1,												// problem status
				((double*) (((int*) dual_result_buffer) + 2)) + 1,								// Obj (UB)
				((double*) (((int*) dual_result_buffer) + 2)) + 2,								// LB
				((double*) (((int*) dual_result_buffer) + 2)) + 3,								// v*
				&global_opt) != 0 ){
				printf("\nProblem while attempting to solve period relaxation of dual scenario subroblem for scenario %s.\n",
					scenarios.buffer + (scenarios.atom_len+1) * *((int*) task_dual_scen_buffer));
				return(-1);
			}
			
			// register worker time
			current_utc_time(&endt);
			*((double*) (((int*) dual_result_buffer) + 2)) = ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
			
			// return results to coordinator
			if( MPI_Send(dual_result_buffer, 1, dualjobresult, 0, RESULT_DUAL_SCEN_INIT, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending initialization result through MPI_Send.\n");
				return(-1);
			}
			
		}else if( task_status.MPI_TAG == DUAL_SCEN_MILP ){		// dual MILP scenario subproblem
			
			// receive scenario and multipliers
			if( MPI_Recv(task_dual_scen_buffer, 1, dualscenjobmsg, 0, task_status.MPI_TAG,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from Coordinator.\n", task_status.MPI_TAG);
				return(-1);
			}
			
			// put scenario into result buffer
			*((int*) dual_result_buffer) = *((int*) task_dual_scen_buffer);
			
			// solve period relaxation with the provided multipliers
			if(global_opt.verbosity >= VERB_WORKER)
				printf("\nP%d: Solving dual scenario subproblem for scenario %s.\n",
					worldrank, scenarios.buffer + (scenarios.atom_len+1) * *((int*) task_dual_scen_buffer));
			if( solve_dual_scenario_subproblem(*(*argv + 1), *(*argv + 3),
				&row_names, &col_names,
				&time_stages, rows_tstage, cols_tstage,
				numfirststagecols, firststagecols, firststagelb, firststageub,
				&scenarios, &scens_parent, scens_branch_tstage,
				scens_start_line, scens_end_line, scens_probability,
				&prob, *((int*) task_dual_scen_buffer), (double*) (((int*) task_dual_scen_buffer) + 1),
				((int*) dual_result_buffer) + 1,												// problem status
				((double*) (((int*) dual_result_buffer) + 2)) + 1,								// Obj (UB)
				((double*) (((int*) dual_result_buffer) + 2)) + 2,								// LB
				((double*) (((int*) dual_result_buffer) + 2)) + 3,								// v*
				&global_opt) != 0 ){
				printf("\nProblem while attempting to solve dual scenario subroblem for scenario %s.\n",
					scenarios.buffer + (scenarios.atom_len+1) * *((int*) task_dual_scen_buffer));
				return(-1);
			}
			
			// register worker time
			current_utc_time(&endt);
			*((double*) (((int*) dual_result_buffer) + 2)) = ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
			
			// return results to coordinator
			if( MPI_Send(dual_result_buffer, 1, dualjobresult, 0, RESULT_DUAL_SCEN_MILP, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending scenario subproblem result through MPI_Send.\n");
				return(-1);
			}
		
		} else if(task_status.MPI_TAG == DUAL_NONSCEN_MILP ){		// dual non-scenario subproblems
			
			// receive LBscen, delayed multipliers, current multipliers and incumbent
			if( MPI_Recv(task_dual_nonscen_buffer, 1, dualnonscenjobmsg, 0, task_status.MPI_TAG,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from Coordinator.\n", task_status.MPI_TAG);
				return(-1);
			}
			
			// put non scenario indicator into result buffer
			*((int*) dual_result_buffer) = NOSCENID;
			
			// solve non scenario dual subproblems
			if(global_opt.verbosity >= VERB_WORKER)
				printf("\nP%d: Solving dual non-scenario subproblems.\n", worldrank);
			if( solve_nonscenario_subproblems(*(*argv + 1),
				&col_names, &row_names, rows_tstage,
				firststageind, numfirststagecols, firststagecols,
				firststagelb, firststageub,
				&prob,
				((double*) (((int*) task_dual_nonscen_buffer) + 1)) + 1,						// delayed multipliers
				((double*) (((int*) task_dual_nonscen_buffer) + 1)) + 1 + numfirststagecols,	// current multipliers
				((double*) (((int*) task_dual_nonscen_buffer) + 1)) + 1 + 2*numfirststagecols,	// incumbent (ucenter)
				((int*) dual_result_buffer) + 1,												// problem status
				((double*) (((int*) dual_result_buffer) + 2)) + 1,								// Obj (UB)
				((double*) (((int*) dual_result_buffer) + 2)) + 2,								// LB
				((double*) (((int*) dual_result_buffer) + 2)) + 3,								// u*
				&global_opt) != 0){
				printf("\nProblem while attempting to solve non-scenario subproblem.\n");
				return(-1);
			}
			
			// add LBscen to non scenario subproblem bounds
			*(((double*) (((int*) dual_result_buffer) + 2)) + 2) +=
				*((double*) (((int*) task_dual_nonscen_buffer) + 1));							// LB for stochastic program
			
			// register worker time
			current_utc_time(&endt);
			*((double*) (((int*) dual_result_buffer) + 2)) = ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
			
			// return results to coordinator
			if( MPI_Send(dual_result_buffer, 1, dualjobresult, 0, RESULT_DUAL_NONSCEN_MILP, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending non-scenario subproblems result through MPI_Send.\n");
				return(-1);
			}
			
		} else if( task_status.MPI_TAG == PRIMAL_FEAS_PROJ ){		// project unfeasible candidate in 1st stage feasible set
			
			// receive infeasible candidate
			if( MPI_Recv(task_primal_buffer, 1, primaljobmsg, 0, task_status.MPI_TAG,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from Coordinator.\n", task_status.MPI_TAG);
				return(-1);
			}
			
			// put candidate ID into result buffer
			*((int*) primal_proj_result_buffer) = *((int*) task_primal_buffer);
			
			// solve projection problem
			if(global_opt.verbosity >= VERB_WORKER)
				printf("\nP%d: Projecting infeasible candidate into 1st stage feasible set.\n", worldrank);
			if( project_infeasible_candidate(*(*argv + 1),
				&col_names, &row_names, rows_tstage,
				firststageind, numfirststagecols, firststagecols,
				firststagelb, firststageub,
				&prob,
				((double*) (((int*) task_primal_buffer) + 2)),								// infeasible candidate
				NULL,																		// problem status
				((double*) (((int*) primal_proj_result_buffer) + 1)) + 1,					// feasible candidate
				&global_opt) != 0){
				printf("\nProblem while attempting to project candidate %d.\n", *((int*) primal_proj_result_buffer));
				return(-1);
			}
			
			// register worker time
			current_utc_time(&endt);
			*((double*) (((int*) primal_proj_result_buffer) + 1)) = ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
			
			// return results to coordinator
			if( MPI_Send(primal_proj_result_buffer, 1, primalprojjobresult, 0,
				RESULT_PRIMAL_FEAS_PROJ, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending primal projection job result through MPI_Send.\n");
				return(-1);
			}
			
		} else if( task_status.MPI_TAG == PRIMAL_SCEN_MILP ){		// primal evaluation of feasible candidate
			
			// receive (candidate, scenario) assignement
			if( MPI_Recv(task_primal_buffer, 1, primaljobmsg, 0, task_status.MPI_TAG,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS ){
				printf("\nError while receiving message %d from Coordinator.\n", task_status.MPI_TAG);
				return(-1);
			}
			
			// put candidate ID and scenario into result buffer
			*((int*) primal_result_buffer) = *((int*) task_primal_buffer);				// candidate ID
			*(((int*) primal_result_buffer) + 1) = *(((int*) task_primal_buffer) + 1);	// scenario number
			
			// evaluating primal candidate
			if(global_opt.verbosity >= VERB_WORKER)
				printf("\nP%d: Evaluating primal candidate for scenario %s.\n",
					worldrank, scenarios.buffer + (scenarios.atom_len+1) * *(((int*) task_primal_buffer) + 1));
			if( evaluate_primal_candidate(*(*argv + 1), *(*argv + 3),
				&row_names, &col_names,
				&time_stages, rows_tstage, cols_tstage,
				numfirststagecols, firststagecols,
				&scenarios, &scens_parent, scens_branch_tstage,
				scens_start_line, scens_end_line, scens_probability,
				&prob,
				*(((int*) task_primal_buffer) + 1),											// scenario number
				(double*) (((int*) task_primal_buffer) + 2),								// candidate
				((int*) primal_result_buffer) + 2,											// problem status
				((double*) (((int*) primal_result_buffer) + 3)) + 1,						// UB
				((double*) (((int*) primal_result_buffer) + 3)) + 2,						// LB
				&global_opt) != 0 ){
				printf("\nProblem while attempting to evaluate primal candidate for scenario %s.\n",
					scenarios.buffer + (scenarios.atom_len+1) * *(((int*) task_primal_buffer) + 1));
				return(-1);
			}
			
			// register worker time
			current_utc_time(&endt);
			*((double*) (((int*) primal_result_buffer) + 3)) = ((double) (endt.tv_sec - startt.tv_sec)) +
				((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC;
			
			// return results to coordinator
			if( MPI_Send(primal_result_buffer, 1, primaljobresult, 0,
				RESULT_PRIMAL_SCEN_MILP, MPI_COMM_WORLD) != MPI_SUCCESS ){
				printf("\nError while sending primal job result through MPI_Send.\n");
				return(-1);
			}
			
		} else if( task_status.MPI_TAG == RELEASE_WORKER ){			// stop worker execution
			
			if(global_opt.verbosity >= VERB_COORDINATOR_ONLY)
				printf("\nWorker %d received termination signal from coordinator.\n", worldrank);
			break;
			
		} else {
			
			printf("\nUnrecognized task received %d while probing for message.\n",
				task_status.MPI_TAG);
			return(-1);
			
		}
		
	}
	
	// Clean heap
	free_options(&global_opt);
	if( free_heap_from_read_SMPS(
		&prob, &row_names, &col_names,
		&time_stages, &rows_tstage, &cols_tstage,
		&scenarios, &scens_parent, &scens_branch_tstage,
		&scens_start_line, &scens_end_line, &scens_probability,
		&firststagecols, &firststagelb, &firststageub, &global_opt) != 0 ) {
		printf("\nError while cleaning the heap from variables storing the SMPS instance. Program will exit.\n");
		return(-1);
	}
	free_message_types(&dualscenjobmsg, &dualnonscenjobmsg, &dualjobresult,
		&primaljobmsg, &primalprojjobresult, &primaljobresult);
	free(task_dual_scen_buffer);
	free(task_dual_nonscen_buffer);
	free(task_primal_buffer);
	free(dual_result_buffer);
	free(primal_proj_result_buffer);
	free(primal_result_buffer);
	
	// Delete the problem and terminate Xpress
	XPRSdestroyprob(prob);
	XPRSfree();
	
	// Return success indicator
	return(0);
}

// Solve a relaxation of scenario subproblem to initialize lower bounds
int solve_init_scenario_subproblem(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, const double *dual_multipliers,
	int *probstatus, double *UB, double *LB, double *firststagesolution, const struct options* global_opt)
{
	// load scenario subproblem
	struct scen_differences restore_point;
	if( load_scenario_subproblem(workdir, rootname, s,
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while loading scenario subproblem for scenario %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// load multipliers
	double *original_obj = (double*) malloc(sizeof(double)*numfirststagecols);
	load_scenario_multipliers(workingproblem, numfirststagecols, firststagecols,
		firststagelb, firststageub, dual_multipliers, original_obj);
	
	// re-enforce cross-scenario bounds on first stage variables (in case they have been changed by the scenario)
	char *bndtype = (char*) malloc(sizeof(char)*numfirststagecols);
	int i;
	for(i = 0; i < numfirststagecols; i++) *(bndtype + i) = 'L';
	XPRSchgbounds(*workingproblem, numfirststagecols, firststagecols, bndtype, firststagelb);
	for(i = 0; i < numfirststagecols; i++) *(bndtype + i) = 'U';
	XPRSchgbounds(*workingproblem, numfirststagecols, firststagecols, bndtype, firststageub);
	free(bndtype);
	
	// allocate relaxed solution vector
	double *relaxed_solution;
	if( (*global_opt).primal_recovery == PRIMAL_RECOVERY_FIFO
		|| (*global_opt).primal_recovery == PRIMAL_RECOVERY_RND
		|| (*global_opt).primal_recovery == PRIMAL_RECOVERY_LIFO ){
		relaxed_solution = (double*) malloc(numfirststagecols*sizeof(double));
	}else if( (*global_opt).primal_recovery == PRIMAL_RECOVERY_IS ){
		relaxed_solution = firststagesolution;
	}else{
		printf("\nUnrecognized primal recovery method %d.\n", (*global_opt).primal_recovery);
		return(-1);
	}
	
	// solve relaxation
	*UB = DECOMP_INFINITY;
	if( (*global_opt).init_type == INIT_TYPE_PERIOD ){
		// read range of delayed rows (they will be ignored from the relaxation)
		int num_delayed_rows, delayed_rows_start;
		char corefilename[200];
		sprintf(corefilename, "%s" SYSTEM_SLASH "%s.mps", workdir, rootname);
		if( read_delayed_rows(corefilename, rownames,
			&num_delayed_rows, &delayed_rows_start) != 0 ){
			printf("\nProblem while reading delayed rows from CORE file:\n\t%s\n", corefilename);
			return(-1);
		}
		//printf("\nNum delayed rows: %d\tStart delayed rows: %d\n", num_delayed_rows, delayed_rows_start);
		// solve period relaxation
		if( solve_sequential_period_relaxation(workingproblem,
			rownames, colnames, num_delayed_rows, delayed_rows_start,
			timestages, rowtstages, coltstages,
			firststageindex, numfirststagecols, firststagecols,
			LB, relaxed_solution, global_opt)  != 0 ){
			printf("\nError while solving period relaxation of scenario subproblem %s.\n",
				(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
			return(-1);
		}
	}else if( (*global_opt).init_type == INIT_TYPE_LINEAR ){
		// solve linear relaxation without delayed constraints and without crossover (if using barrier)
		int crossover_setting;
		char mipoptimizeflags[3];
		XPRSgetintcontrol(*workingproblem, XPRS_CROSSOVER, &crossover_setting);
		XPRSsetintcontrol(*workingproblem, XPRS_CROSSOVER, 0);
		strcpy(mipoptimizeflags, (*global_opt).scenprob_options.lp_method);
		strcat(mipoptimizeflags, "l");
		if( XPRSmipoptimize(*workingproblem, mipoptimizeflags) != 0 ){
			printf("Xpress failed at solving initialization for scenario %s, non-zero returned from XPRSmipoptimize.",
				(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
			return(-1);
		}
		XPRSgetintattrib(*workingproblem, XPRS_MIPSTATUS, probstatus);
		XPRSpostsolve(*workingproblem);
		if( !( *probstatus == XPRS_MIP_OPTIMAL
			|| *probstatus == XPRS_MIP_SOLUTION
			|| *probstatus == XPRS_MIP_LP_OPTIMAL ) ){
			printf("Xpress failed at solving scenario subproblem for scenario %s."
				"\nProblem status (MIPSTATUS):\t", (*scenarios).buffer + ((*scenarios).atom_len+1)*s);
			switch(*probstatus){
				case XPRS_MIP_NOT_LOADED:
					printf("XPRS_MIP_NOT_LOADED");
					break;
				case XPRS_MIP_LP_NOT_OPTIMAL:
					printf("XPRS_MIP_LP_NOT_OPTIMAL");
					break;
				case XPRS_MIP_NO_SOL_FOUND:
					printf("XPRS_MIP_NO_SOL_FOUND");
					break;
				case XPRS_MIP_INFEAS:
					printf("XPRS_MIP_INFEAS");
					break;
				case XPRS_MIP_UNBOUNDED:
					printf("XPRS_MIP_UNBOUNDED");
					break;
				default: printf("Unknown status");
			}
			export_troublemaker(workdir, "unsolved_initialization_subprob.lp", "l", workingproblem);
			return(-1);
		}
		// parse results to output buffers
		double *lp_full_sol;
		XPRSgetdblattrib(*workingproblem, XPRS_BESTBOUND, LB);
		lp_full_sol = (double*) malloc(sizeof(double)*(*colnames).num_elems);
		XPRSgetlpsol(*workingproblem, lp_full_sol, NULL, NULL, NULL);
		int j;
		for(j = 0; j < numfirststagecols; j++){
			*(relaxed_solution + j) = *(lp_full_sol + *(firststagecols + j));
		}
		free(lp_full_sol);
		// post solve the matrix for the next subproblem
		XPRSpostsolve(*workingproblem);
	}else{
		printf("\nUnrecognized initialization method.\n");
		return(-1);
	}
	
	// project relaxed solution and free malloc'd space if necessary
	if( (*global_opt).primal_recovery == PRIMAL_RECOVERY_FIFO
		|| (*global_opt).primal_recovery == PRIMAL_RECOVERY_RND
		|| (*global_opt).primal_recovery == PRIMAL_RECOVERY_LIFO ){
		if( project_infeasible_candidate(workdir, colnames, rownames, rowtstages,
			firststageindex, numfirststagecols, firststagecols,
			firststagelb, firststageub,
			workingproblem, relaxed_solution,
			NULL, firststagesolution, global_opt) != 0){
			printf("\nError while projecting initialization relaxation solution of scenario subproblem %s"
				"\nonto the first stage feasible set.\n",
				(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
			return(-1);
		}
		free(relaxed_solution);
	}
	
	// remove multipliers
	unload_scenario_multipliers(workingproblem, numfirststagecols, firststagecols, original_obj);
	free(original_obj);
	
	// restore worker subproblem
	if( restore_core_problem(
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while restoring CORE problem from scenario subproblem %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	//char fname[100];
	//strcpy(fname, "RelaxedSolution");
	//strcat(fname, (*scenarios).buffer + ((*scenarios).atom_len + 1)*s);
	//strcat(fname, ".csv");
	//FILE *relaxsol = fopen(fname, "w");
	//int i;
	//for(i = 0; i < numfirststagecols; i++){
	//	fprintf(relaxsol, "%f\n", *(relaxed_solution + i));
	//}
	//fclose(relaxsol);
	
	// Return success indicator
	return(0);
}

int solve_dual_scenario_subproblem(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	int numfirststagecols, const int *firststagecols, const double *firststagelb, const double *firststageub,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, const double *dual_multipliers,
	int *probstatus, double *UB, double *LB, double *firststagesolution, const struct options* global_opt)
{
	/*
	// write received multipliers
	{
		int j;
		char multiplierfile[100];
		FILE *out;
		sprintf(multiplierfile, "xlast_%d.csv", s);
		out = fopen(multiplierfile, "w");
		fprintf(out, "j,x\n");
		for(j = 0; j < numfirststagecols; j++){
			fprintf(out, "%d,%f\n", j, *(dual_multipliers + j));
		}
		fclose(out);
	}
	*/
	
	// load scenario subproblem
	struct scen_differences restore_point;
	if( load_scenario_subproblem(workdir, rootname, s,
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while loading scenario subproblem for scenario %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// load multipliers
	double *original_obj = (double*) malloc(sizeof(double)*numfirststagecols);
	load_scenario_multipliers(workingproblem, numfirststagecols, firststagecols,
		firststagelb, firststageub, dual_multipliers, original_obj);
	
	// re-enforce cross-scenario bounds on first stage variables (in case they have been changed by the scenario)
	char *bndtype = (char*) malloc(sizeof(char)*numfirststagecols);
	int i;
	for(i = 0; i < numfirststagecols; i++) *(bndtype + i) = 'L';
	XPRSchgbounds(*workingproblem, numfirststagecols, firststagecols, bndtype, firststagelb);
	for(i = 0; i < numfirststagecols; i++) *(bndtype + i) = 'U';
	XPRSchgbounds(*workingproblem, numfirststagecols, firststagecols, bndtype, firststageub);
	free(bndtype);
	
	// solve subproblem
	if( XPRSmipoptimize(*workingproblem, (*global_opt).scenprob_options.lp_method) != 0 ){
		printf("Xpress failed at solving scenario subproblem for scenario %s, non-zero returned from XPRSmipoptimize.",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	//XPRSgetintattrib(workingproblem, XPRS_MIPSTATUS, (int*) result_buffer + 1);
	XPRSgetintattrib(*workingproblem, XPRS_MIPSTATUS, probstatus);
	
	if( !( *probstatus == XPRS_MIP_OPTIMAL || *probstatus == XPRS_MIP_SOLUTION ) ){
		printf("Xpress failed at solving scenario subproblem for scenario %s."
			"\nProblem status (MIPSTATUS):\t", (*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		switch(*probstatus){
			case XPRS_MIP_NOT_LOADED:
				printf("XPRS_MIP_NOT_LOADED");
				break;
			case XPRS_MIP_LP_NOT_OPTIMAL:
				printf("XPRS_MIP_LP_NOT_OPTIMAL");
				break;
			case XPRS_MIP_LP_OPTIMAL:
				printf("XPRS_MIP_LP_OPTIMAL");
				break;
			case XPRS_MIP_NO_SOL_FOUND:
				printf("XPRS_MIP_NO_SOL_FOUND");
				break;
			case XPRS_MIP_INFEAS:
				printf("XPRS_MIP_INFEAS");
				break;
			case XPRS_MIP_UNBOUNDED:
				printf("XPRS_MIP_UNBOUNDED");
				break;
			default: printf("Unknown status");
		}
		export_troublemaker(workdir, "unsolved_dual_scenario_subprob.lp", "l", workingproblem);
		return(-1);
	}
	
	// postsolve if necessary
	if( XPRSpostsolve(*workingproblem) != 0 ){
		printf("Xpress failed at postsolving scenario subproblem for scenario %s.",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// parse results to output buffers
	double *mip_full_sol;
	//XPRSgetdblattrib(prob, XPRS_MIPOBJVAL, (double*) ((int*) result_buffer + 2));
	XPRSgetdblattrib(*workingproblem, XPRS_MIPOBJVAL, UB);
	//XPRSgetdblattrib(prob, XPRS_BESTBOUND, (double*) ((int*) result_buffer + 2) + 1);
	XPRSgetdblattrib(*workingproblem, XPRS_BESTBOUND, LB);
	mip_full_sol = (double*) malloc(sizeof(double)*(*colnames).num_elems);
	XPRSgetmipsol(*workingproblem, mip_full_sol, NULL);
	int j;
	for(j = 0; j < numfirststagecols; j++){
		//*(u + i) = *(mip_full_sol + *(firststagecols + i));
		*(firststagesolution + j) = *(mip_full_sol + *(firststagecols + j));
	}
	free(mip_full_sol);
	
	// remove multipliers
	unload_scenario_multipliers(workingproblem, numfirststagecols, firststagecols, original_obj);
	free(original_obj);
	
	// restore worker subproblem
	if( restore_core_problem(
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while restoring CORE problem from scenario subproblem %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

int solve_nonscenario_subproblems(const char *workdir, const struct string_buffer *colnames,
	const struct string_buffer *rownames, const int *rowtstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	XPRSprob *workingproblem, const double *msum_delayed_multipliers,
	const double *msum_current_multipliers, const double *ucenter,
	int *probstatus, double *UB, double *LB, double *solution, const struct options *global_opt)
{
	/*
	// write multipliers to files
	int i;
	FILE *nonscenfile = fopen("NonScenAssignement.csv", "w");
	fprintf(nonscenfile, "COLUMN,msumx_delayed,msumx_current,ucenter\n");
	for(i = 0; i < numfirststagecols; i++){
		fprintf(nonscenfile, "%s,%f,%f,%f\n",
			(*colnames).buffer + ( (*colnames).atom_len + 1 ) * *(firststagecols + i),
			*(msum_delayed_multipliers + i), *(msum_current_multipliers + i),
			*(ucenter + i));
	}
	fclose(nonscenfile);
	*/
	
	// create nonscenario problem
	XPRSprob f0prob;
	XPRScreateprob(&f0prob);
	XPRSsetintcontrol(f0prob, XPRS_SCALING, 0);				// avoids automatic scaling when reading the problem
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		XPRSaddcbmessage(f0prob, Message, NULL, 0);			// defines callback function for printing log from Xpress (usefull only in Windows)
	#endif
	
	// load xpress parameters
	{
		int i;
		for(i = 0; i < (*global_opt).periodprob_options.num_intoptions; i++){
			XPRSsetintcontrol(f0prob, *((*global_opt).periodprob_options.intoptions + i),
				*((*global_opt).periodprob_options.intoptions_value + i));
		}
		for(i = 0; i < (*global_opt).periodprob_options.num_dbloptions; i++){
			XPRSsetdblcontrol(f0prob, *((*global_opt).periodprob_options.dbloptions + i),
				*((*global_opt).periodprob_options.dbloptions_value + i));
		}
	}
	
	// override verbosity level
	if((*global_opt).verbosity > VERB_WORKER){
		XPRSsetintcontrol(f0prob, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(f0prob, XPRS_OUTPUTLOG, 4);
	}
	
	// formulate first stage problem
	load_first_stage_matrix(workingproblem, rownames, rowtstages, firststageindex,
		numfirststagecols, firststagecols, &f0prob, "f0_problem");
	
	// load delayed multipliers
	load_f0_lin_objective(numfirststagecols, firststagelb, firststageub,
		msum_delayed_multipliers, &f0prob);
	
	// solve subproblem
	if( XPRSmipoptimize(f0prob, (*global_opt).periodprob_options.lp_method) != 0 ){
		printf("\nXpress failed at solving non-scenario subproblem with delayed multipliers,"
			"\nnon-zero returned from XPRSmipoptimize.");
		return(-1);
	}
	XPRSgetintattrib(f0prob, XPRS_MIPSTATUS, probstatus);
	if( !(*probstatus == XPRS_MIP_OPTIMAL || *probstatus == XPRS_MIP_SOLUTION) ){
		printf("\nXpress failed at solving non-scenario subproblem with delayed multipliers."
			"\nProblem status (MIPSTATUS):\t");
		switch(*probstatus){
			case XPRS_MIP_NOT_LOADED:
				printf("XPRS_MIP_NOT_LOADED");
				break;
			case XPRS_MIP_LP_NOT_OPTIMAL:
				printf("XPRS_MIP_LP_NOT_OPTIMAL");
				break;
			case XPRS_MIP_LP_OPTIMAL:
				printf("XPRS_MIP_LP_OPTIMAL");
				break;
			case XPRS_MIP_NO_SOL_FOUND:
				printf("XPRS_MIP_NO_SOL_FOUND");
				break;
			//case XPRS_MIP_SOLUTION:
			//	printf("XPRS_MIP_SOLUTION");
			//	break;
			case XPRS_MIP_INFEAS:
				printf("XPRS_MIP_INFEAS");
				break;
			case XPRS_MIP_UNBOUNDED:
				printf("XPRS_MIP_UNBOUNDED");
				break;
			default: printf("Unknown status");
		}
		export_troublemaker(workdir, "unsolved_dual_nonscenario_linsubprob.lp", "l", &f0prob);
		return(-1);
	}
	
	// parse objective into result_buffer
	XPRSgetdblattrib(f0prob, XPRS_MIPOBJVAL, UB);
	XPRSgetdblattrib(f0prob, XPRS_BESTBOUND, LB);
	
	// what to next do depends on the smoothing parameter
	if( fabs((*global_opt).dual_f0_mu) < DBL_EPSILON ){
		
		// load current multiplers
		clean_f0_objective(&f0prob);
		load_f0_lin_objective(numfirststagecols, firststagelb, firststageub,
			msum_current_multipliers, &f0prob);
		
		// resolve subproblem
		if( XPRSmipoptimize(f0prob, (*global_opt).periodprob_options.lp_method) != 0 ){
			printf("\nXpress failed at solving non-scenario subproblem with current multipliers,"
				"\nnon-zero returned from XPRSmipoptimize.");
			return(-1);
		}
		XPRSgetintattrib(f0prob, XPRS_MIPSTATUS, probstatus);
		if( !(*probstatus == XPRS_MIP_OPTIMAL || *probstatus == XPRS_MIP_SOLUTION) ){
			printf("\nXpress failed at solving non-scenario subproblem with current multipliers."
				"\nProblem status (MIPSTATUS):\t");
			switch(*probstatus){
				case XPRS_MIP_NOT_LOADED:
					printf("XPRS_MIP_NOT_LOADED");
					break;
				case XPRS_MIP_LP_NOT_OPTIMAL:
					printf("XPRS_MIP_LP_NOT_OPTIMAL");
					break;
				case XPRS_MIP_LP_OPTIMAL:
					printf("XPRS_MIP_LP_OPTIMAL");
					break;
				case XPRS_MIP_NO_SOL_FOUND:
					printf("XPRS_MIP_NO_SOL_FOUND");
					break;
				//case XPRS_MIP_SOLUTION:
				//	printf("XPRS_MIP_SOLUTION");
				//	break;
				case XPRS_MIP_INFEAS:
					printf("XPRS_MIP_INFEAS");
					break;
				case XPRS_MIP_UNBOUNDED:
					printf("XPRS_MIP_UNBOUNDED");
					break;
				default: printf("Unknown status");
			}
			export_troublemaker(workdir, "unsolved_dual_nonscenario_linsubprob.lp", "l", &f0prob);
			return(-1);
		}
		
		// parse subgradient into result_buffer
		XPRSgetmipsol(f0prob, solution, NULL);
		
	}else if( (*global_opt).dual_f0_mu > DBL_EPSILON ){
		
		// load center
		double quadoffset;
		clean_f0_objective(&f0prob);
		load_f0_quad_objective(0, colnames,
			numfirststagecols, firststagecols, firststagelb, firststageub,
			msum_current_multipliers, (*global_opt).dual_f0_mu, ucenter,
			&f0prob, &quadoffset);
		
		// solve subproblem
		if( XPRSlpoptimize(f0prob, "b") != 0 ){
			printf("\nXpress failed at solving non-scenario quadratic subproblem, non-zero returned from XPRSlpoptimize.");
			return(-1);
		}
		XPRSgetintattrib(f0prob, XPRS_LPSTATUS, probstatus);
		if(*probstatus != XPRS_LP_OPTIMAL){
			printf("\nXpress failed at solving non-scenario quadratic subproblem."
				"\nProblem status (LPSTATUS):\t");
			switch(*probstatus){
				case XPRS_LP_UNSTARTED:
					printf("XPRS_LP_UNSTARTED");
					break;
				case XPRS_LP_INFEAS:
					printf("XPRS_LP_INFEAS");
					break;
				case XPRS_LP_CUTOFF:
					printf("XPRS_LP_CUTOFF");
					break;
				case XPRS_LP_UNFINISHED:
					printf("XPRS_LP_UNFINISHED");
					break;
				case XPRS_LP_UNBOUNDED:
					printf("XPRS_LP_UNBOUNDED");
					break;
				case XPRS_LP_CUTOFF_IN_DUAL:
					printf("XPRS_LP_CUTOFF_IN_DUAL");
					break;
				case XPRS_LP_UNSOLVED:
					printf("XPRS_LP_UNSOLVED");
					break;
				case XPRS_LP_NONCONVEX:
					printf("XPRS_LP_NONCONVEX");
					break;
				default: printf("Unknown status");
			}
			export_troublemaker(workdir, "unsolved_dual_nonscenario_quadsubprob.lp", "l", &f0prob);
			return(-1);
		}
		
		// parse gradient into result_buffer
		XPRSgetlpsol(f0prob, solution, NULL, NULL, NULL);
		
	}else{
		
		printf("\nNegative value for the smoothing parameter of f0_mu: %f\n", (*global_opt).dual_f0_mu);
		return(-1);
		
	}
	
	// destroy first stage problem
	XPRSdestroyprob(f0prob);
	
	// Return success indicator
	return(0);
}

int project_infeasible_candidate(const char *workdir, const struct string_buffer *colnames,
	const struct string_buffer *rownames, const int *rowtstages,
	const int firststageindex, const int numfirststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub,
	XPRSprob *workingproblem, const double *infeasible_candidate,
	int *probstatus, double *solution, const struct options* global_opt)
{
	// Iterators
	int i;
	
	// create nonscenario problem
	XPRSprob f0prob;
	XPRScreateprob(&f0prob);
	XPRSsetintcontrol(f0prob, XPRS_SCALING, 0);				// avoids automatic scaling when reading the problem
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		XPRSaddcbmessage(f0prob, Message, NULL, 0);			// defines callback function for printing log from Xpress (usefull only in Windows)
	#endif
	
	// load xpress parameters
	for(i = 0; i < (*global_opt).periodprob_options.num_intoptions; i++){
		XPRSsetintcontrol(f0prob, *((*global_opt).periodprob_options.intoptions + i),
			*((*global_opt).periodprob_options.intoptions_value + i));
	}
	for(i = 0; i < (*global_opt).periodprob_options.num_dbloptions; i++){
		XPRSsetdblcontrol(f0prob, *((*global_opt).periodprob_options.dbloptions + i),
			*((*global_opt).periodprob_options.dbloptions_value + i));
	}
	
	// override verbosity level
	if((*global_opt).verbosity > VERB_WORKER){
		XPRSsetintcontrol(f0prob, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(f0prob, XPRS_OUTPUTLOG, 4);
	}
	
	// override tolerances to 2 orders of magnited below tolerance for scenario subproblems
	double scenprob_FEASTOL, scenprob_MIPTOL;
	/*
	for(i = 0; i < (*global_opt).scenprob_options.num_dbloptions; i++){
		if( *((*global_opt).scenprob_options.dbloptions + i) == XPRS_MIPTOL ){
			break;
		}
	}
	if( i < (*global_opt).scenprob_options.num_dbloptions ){
		scenprob_MIPTOL = *((*global_opt).scenprob_options.dbloptions_value + i);
	}else{ // MIPTOL option not specified for scenario subproblems -> reducing default value
		scenprob_MIPTOL = uXPRS_MIPTOL_DEFAULT;
	}
	*/
	XPRSgetdblcontrol(*workingproblem, XPRS_FEASTOL, &scenprob_FEASTOL);
	XPRSgetdblcontrol(*workingproblem, XPRS_MIPTOL, &scenprob_MIPTOL);
	XPRSsetdblcontrol(f0prob, XPRS_FEASTOL, 0.01 * scenprob_FEASTOL);
	XPRSsetdblcontrol(f0prob, XPRS_MIPTOL, 0.01 * scenprob_MIPTOL);
	
	// formulate first stage problem
	load_first_stage_matrix(workingproblem, rownames, rowtstages, firststageindex,
		numfirststagecols, firststagecols, &f0prob, "f0_problem");
	
	// get global variables (binaries)
	int nglents, nsets;
	char *qgtype = (char*) malloc(numfirststagecols*sizeof(char));
	int *mgcols = (int*) malloc(numfirststagecols*sizeof(int));
	XPRSgetglobal(f0prob, &nglents, &nsets, qgtype, mgcols,
		NULL, NULL, NULL, NULL, NULL);
	
	// adjust infeasible candidate according to the specified midpoint for rounding binaries
	double *infeasible_candidate_adjusted = (double*) malloc(numfirststagecols*sizeof(double));
	memcpy(infeasible_candidate_adjusted, infeasible_candidate, numfirststagecols*sizeof(double));
														// here I created an internal copy to avoid modifying the extenal candidate
	if( fabs((*global_opt).proj_round_midpoint - .5) > DBL_EPSILON ){
		for(i = 0; i < nglents; i++){
			if( *(qgtype + i) == 'B' ){
				if( *(infeasible_candidate_adjusted + *(mgcols + i)) <= (*global_opt).proj_round_midpoint ){
					*(infeasible_candidate_adjusted + *(mgcols + i)) *=
						1.0 + (0.5 - (*global_opt).proj_round_midpoint)/(*global_opt).proj_round_midpoint;
				} else {
					*(infeasible_candidate_adjusted + *(mgcols + i)) =
						(1.0 - (0.5 - (*global_opt).proj_round_midpoint)/(1-(*global_opt).proj_round_midpoint)) *
						*(infeasible_candidate_adjusted + *(mgcols + i)) +
						(0.5 - (*global_opt).proj_round_midpoint)/(1-(*global_opt).proj_round_midpoint);
				}
			}
		}
	}
	
	// load infeasible candidate in the objective
	double quadoffset;
	double *zero_multipliers = (double*) calloc(numfirststagecols, sizeof(double));
	clean_f0_objective(&f0prob);
	load_f0_quad_objective(1, colnames,
		numfirststagecols, firststagecols, firststagelb, firststageub,
		zero_multipliers, 1, infeasible_candidate_adjusted, &f0prob, &quadoffset);
	free(zero_multipliers);
	
	// export problem (debugging)
	//export_troublemaker(workdir, "unsolved_feasibility_projection_subprob.lp", "l", &f0prob);
	
	// solve subproblem
	XPRSmipoptimize(f0prob, "b");
	//if( XPRSmipoptimize(f0prob, "b") != 0 ){
	//	printf("\nXpress failed at solving feasibility recovery quadratic subproblem, non-zero returned from XPRSmipoptimize.");
	//	export_troublemaker(workdir, "unsolved_projection_subprob.lp", "l", &f0prob);
	//	return(-1);
	//}
	int localprobstatus;
	XPRSgetintattrib(f0prob, XPRS_MIPSTATUS, &localprobstatus);
	if( probstatus != NULL) *probstatus = localprobstatus;
	if(!(localprobstatus == XPRS_MIP_OPTIMAL || localprobstatus == XPRS_MIP_SOLUTION)){
		printf("\nXpress failed at solving feasibility recovery subproblem."
			"\nProblem status (MIPSTATUS):\t");
		switch(localprobstatus){
			case XPRS_MIP_NOT_LOADED:
				printf("XPRS_MIP_NOT_LOADED");
				break;
			case XPRS_MIP_LP_NOT_OPTIMAL:
				printf("XPRS_MIP_LP_NOT_OPTIMAL");
				break;
			case XPRS_MIP_LP_OPTIMAL:
				printf("XPRS_MIP_LP_OPTIMAL");
				break;
			case XPRS_MIP_NO_SOL_FOUND:
				printf("XPRS_MIP_NO_SOL_FOUND");
				break;
			//case XPRS_MIP_SOLUTION:
			//	printf("XPRS_MIP_SOLUTION");
			//	break;
			case XPRS_MIP_INFEAS:
				printf("XPRS_MIP_INFEAS");
				break;
			case XPRS_MIP_UNBOUNDED:
				printf("XPRS_MIP_UNBOUNDED");
				break;
			default: printf("Unknown status %d", localprobstatus);
		}
		export_troublemaker(workdir, "unsolved_feasibility_projection_subprob.lp", "l", &f0prob);
		return(-1);
	}
	
	// debugging
	//return(-1);
	
	// parse solution into result_buffer
	XPRSgetmipsol(f0prob, solution, NULL);
	
	// round binaries and integers
	for(i = 0; i < nglents; i++){
		if( *(qgtype + i) == 'B' || *(qgtype + i) == 'I' ){
			*(solution + *(mgcols + i)) = round(*(solution + *(mgcols + i)));
		}
	}
	
	// resolve problem with rounded integers if necessary
	if(numfirststagecols > nglents){
		
		// fixing all globals to their current values
		char *qbtype = (char*) malloc(sizeof(char)*nglents);
		double *bnd = (double*) malloc(sizeof(double)*nglents);
		for(i = 0; i < nglents; i++) *(qbtype + i) = 'B';
		for(i = 0; i < nglents; i++) *(bnd + i) = *(solution + *(mgcols + i));
		XPRSchgbounds(f0prob, nglents, mgcols, qbtype, bnd);
		free(qbtype);
		free(bnd);
		
		// resolve program with continuous variables only
		XPRSmipoptimize(f0prob, "b");
		XPRSgetintattrib(f0prob, XPRS_MIPSTATUS, &localprobstatus);
		if( probstatus != NULL) *probstatus = localprobstatus;
		if(!(localprobstatus == XPRS_MIP_OPTIMAL || localprobstatus == XPRS_MIP_SOLUTION)){
			printf("\nXpress failed at solving feasibility recovery subproblem."
				"\nProblem status (MIPSTATUS):\t");
			switch(localprobstatus){
				case XPRS_MIP_NOT_LOADED:
					printf("XPRS_MIP_NOT_LOADED");
					break;
				case XPRS_MIP_LP_NOT_OPTIMAL:
					printf("XPRS_MIP_LP_NOT_OPTIMAL");
					break;
				case XPRS_MIP_LP_OPTIMAL:
					printf("XPRS_MIP_LP_OPTIMAL");
					break;
				case XPRS_MIP_NO_SOL_FOUND:
					printf("XPRS_MIP_NO_SOL_FOUND");
					break;
				//case XPRS_MIP_SOLUTION:
				//	printf("XPRS_MIP_SOLUTION");
				//	break;
				case XPRS_MIP_INFEAS:
					printf("XPRS_MIP_INFEAS");
					break;
				case XPRS_MIP_UNBOUNDED:
					printf("XPRS_MIP_UNBOUNDED");
					break;
				default: printf("Unknown status %d", localprobstatus);
			}
			export_troublemaker(workdir, "unsolved_feasibility_second_projection_subprob.lp", "l", &f0prob);
			return(-1);
		}
		
		// store new solution
		XPRSgetmipsol(f0prob, solution, NULL);
		
	}
	
	// free memory
	free(qgtype);
	free(mgcols);
	free(infeasible_candidate_adjusted);
	
	// destroy first stage problem
	XPRSdestroyprob(f0prob);
	
	// Return success indicator
	return(0);
}

int evaluate_primal_candidate(const char *workdir, const char *rootname,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int numfirststagecols, const int *firststagecols,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *workingproblem, const int s, double *first_stage_candidate,
	int *probstatus, double *UB, double *LB, const struct options* global_opt)
{
	// load scenario subproblem
	struct scen_differences restore_point;
	if( load_scenario_subproblem(workdir, rootname, s,
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while loading scenario subproblem for scenario %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// fix first stage variables
	double *firststageLB = malloc(sizeof(double)*numfirststagecols);
	double *firststageUB = malloc(sizeof(double)*numfirststagecols);
	fix_first_stage_vars(workingproblem, numfirststagecols, firststagecols,
		first_stage_candidate, firststageLB, firststageUB);
	
	// solve subproblem
	if( XPRSmipoptimize(*workingproblem, (*global_opt).scenprob_options.lp_method) != 0 ){
		printf("Xpress failed at solving candidate evaluation for scenario %s, non-zero returned from XPRSmipoptimize.",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	XPRSgetintattrib(*workingproblem, XPRS_MIPSTATUS, probstatus);
	
	if( !( *probstatus == XPRS_MIP_OPTIMAL || *probstatus == XPRS_MIP_SOLUTION ) ){
		printf("Xpress failed at solving candidate evaluation for scenario %s."
			"\nProblem status (MIPSTATUS):\t", (*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		switch(*probstatus){
			case XPRS_MIP_NOT_LOADED:
				printf("XPRS_MIP_NOT_LOADED");
				break;
			case XPRS_MIP_LP_NOT_OPTIMAL:
				printf("XPRS_MIP_LP_NOT_OPTIMAL");
				break;
			case XPRS_MIP_LP_OPTIMAL:
				printf("XPRS_MIP_LP_OPTIMAL");
				break;
			case XPRS_MIP_NO_SOL_FOUND:
				printf("XPRS_MIP_NO_SOL_FOUND");
				break;
			case XPRS_MIP_INFEAS:
				printf("XPRS_MIP_INFEAS");
				break;
			case XPRS_MIP_UNBOUNDED:
				printf("XPRS_MIP_UNBOUNDED");
				break;
			default: printf("Unknown status");
		}
		export_troublemaker(workdir, "unsolved_primal_evaluation_subprob.lp", "l", workingproblem);
		return(-1);
	}
	
	// postsolve if necessary
	if( XPRSpostsolve(*workingproblem) != 0 ){
		printf("Xpress failed at postsolving candidate evaluation for scenario %s.",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// parse results to output buffers
	XPRSgetdblattrib(*workingproblem, XPRS_MIPOBJVAL, UB);
	XPRSgetdblattrib(*workingproblem, XPRS_BESTBOUND, LB);
	
	// unfix first stage variables
	unfix_first_stage_vars(workingproblem, numfirststagecols, firststagecols, firststageLB, firststageUB);
	free(firststageLB);
	free(firststageUB);
	
	// restore worker subproblem
	if( restore_core_problem(
		rownames, colnames,
		timestages, rowtstages, coltstages,
		scenarios, scenparent, scenbranchtstage,
		scenstartline, scenendline, scenprobability,
		workingproblem, &restore_point, global_opt) != 0){
		printf("\nError while restoring CORE problem from scenario subproblem %s.\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*s);
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

void export_troublemaker(const char *workdir, const char *filename,
	const char *format, XPRSprob *problem)
{
	// get rank
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	// collect the parts of the full name
	char writeto[1000];
	sprintf(writeto, "%s" SYSTEM_SLASH "P%d_%s", workdir, world_rank, filename);
	
	// issue write command
	XPRSwriteprob(*problem, writeto, format);
	
	// print confirmation
	printf("\nUnsolved problem written to:\n\t%s\n", writeto);
}