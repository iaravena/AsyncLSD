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

#ifndef COORDINATOR_INCLUDED
	
	// define COORDINATOR_INCLUDED to avoid loading the header again
	#define COORDINATOR_INCLUDED
	
	// Constants
	#define RANDOM_SCENARIO		-1
	
	// Metaprocesses
	extern int allocate_memory_space_async_alg(const int numscenarios, const int numfirststagecols,
		const int initial_numcandidates, const int primal_queue_size, double **evaluation_timestamp,
		void **x_evaluated, double ***x_evaluated_scen_dbl, void **v_evaluated, void ***v_evaluated_scen,
		double ***v_evaluated_scen_dbl, void **u_evaluated, double **u_evaluated_dbl, void **x_current,
		void ***x_current_scen, double ***x_current_scen_dbl, void **x_nonscen, double **LBscen,
		double **msumx_evaluated, double **msumx_current, double **ucenter,
		void ***candidates, int **candidate_status, double **candidate_UB, int **candidate_numevals,
		double ***candidate_UB_scen, int **primal_queue_candidate, int **primal_queue_scenario);
	extern void free_memory_space_async_alg(double **evaluation_timestamp,
		void **x_evaluated, double ***x_evaluated_scen_dbl, void **v_evaluated, void ***v_evaluated_scen,
		double ***v_evaluated_scen_dbl, void **u_evaluated, void **x_current,
		void ***x_current_scen, double ***x_current_scen_dbl, void **x_nonscen,
		void ***candidates, int **candidate_status, double **candidate_UB, int **candidate_numevals,
		double ***candidate_UB_scen, int **primal_queue_candidate, int **primal_queue_scenario);
	extern int initialize_subgradients_and_bounds(
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
		struct timespec *startt, FILE *iterfile, const struct options *opt);
	
	// Dual jobs
	extern int update_and_launch_dual_scen_subprob(const int k, const int scenario,
		const double gnorm2, const double LB, const double UB,
		const int numscenarios, const double *scenprobabilities,							// required for stepsize computation
		const int numfirststagecols, const double *update_scaling,
		void **x_current_scen, double **x_current_scen_dbl,
		double *msumx_current, double **v_evaluated_scen_dbl, double *u_evaluated_dbl,		// required for updating
		const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *dualscenjobmsg, struct timespec *startt, const struct options *opt);
	extern int posprocess_dual_scen_subprob(const int k, const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, double *evaluation_timestamp, double **x_evaluated_scen_dbl,
		double *LBscen, double *msumx_evaluated, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
		double *ucenter, MPI_Status *result_status, struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *dualjobresult, int *numcandidates, int *maxnumcandidates, void ***candidates,
		int **candidate_status, double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
		struct timespec *startt, FILE *iterfile, const struct options *opt);
	extern int launch_dual_nonscen_subprob(const int numfirststagecols, const void *x_nonscen,
		const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *dualnonscenjobmsg, struct timespec *startt, const struct options *opt);
	extern int posprocess_dual_nonscen_subprob(const int k, const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, const double *update_scaling,
		void *u_evaluated, double *u_evaluated_dbl, double **v_evaluated_scen_dbl,
		const double UBk, double *LBaux, double *LBk, double *gnorm2, MPI_Status *result_status,
		struct slave *workers, MPI_Request *work_requests, MPI_Datatype *dualjobresult,
		struct timespec *startt, FILE *iterfile, const struct options *opt);
	
	// Primal jobs
	extern int launch_primal_candidate_evaluation(
		const int numfirststagecols, const int numcandidates, void **candidates,
		const int primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
		const int *primal_queue_candidate, const int *primal_queue_scenario,
		const int assigned_worker, int *previous_primal_worker,
		struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *primaljobmsg, struct timespec *startt, const struct options *opt);
	extern int posprocess_primal_candidate_evaluation(const int k,
		const int numscenarios, const double *scenprobabilities,
		void **candidates, int *candidate_status, int *candidate_numevals,
		double *candidate_UB, double **candidate_UB_scen,
		int *best_candidate, double *UBk, const double LBk,
		const int numfirststagecols, double *ucenter,
		MPI_Status *result_status, struct slave *workers,
		MPI_Request *work_requests, MPI_Datatype *primaljobresult,
		struct timespec *startt, FILE *iterfile, FILE *candidatefile, const struct options *opt);
	extern int generate_and_project_IS_primal_candidate(
		const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
		const int assigned_worker, struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *primaljobmsg, struct timespec *startt, const struct options *opt);
	extern int posprocess_primal_feasibility_projection(const int k, const int critical_queue_lenght,
		const int numscenarios, const int numfirststagecols,
		int *numcandidates, int *maxnumcandidates, void ***candidates, int **candidate_status,
		double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
		int *primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
		int *primal_queue_last, int **primal_queue_candidate, int **primal_queue_scenario,
		int *num_projections, double *av_time_projection, double *stddev_time_projection,
		MPI_Status *result_status, struct slave *workers, MPI_Request *work_requests,
		MPI_Datatype *primalprojjobresult, struct timespec *startt, FILE *iterfile,
		const struct options *opt);
	extern int select_next_primal_candidate(const int numcandidates, int *candidate_status,
		int *next_candidate_ID, const struct options *opt);
	extern void generate_new_IS_primal_candidate(const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, void **v_evaluated_scen, double **v_evaluated_scen_dbl,
		double *new_candidate, const struct options *opt);
	extern int add_primal_candidate(const int numscenarios,
		const int numfirststagecols, const double *new_candidate,
		int *numcandidates, int *maxnumcandidates, void ***candidates, int **candidate_status,
		double **candidate_UB, int **candidate_numevals, double ***candidate_UB_scen,
		const struct options *opt);
	extern int activate_primal_candidate(const int candidate_ID, const int numscenarios,
		const int numcandidates, int *candidate_status, double **candidate_UB_scen,
		int *primal_queue_size, int *primal_queue_numjobs, int *primal_queue_first,
		int *primal_queue_last, int **primal_queue_candidate, int **primal_queue_scenario);
	extern void consolidate_primal_candidate(const int candidate_ID,
		const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, void **candidates, int *candidate_status,
		double *candidate_UB, double **candidate_UB_scen, FILE *candidatesfile);
	extern void free_memory_candidates(const int numcandidates, int *candidate_status,
		void **candidates, double **candidate_UB_scen);
	
	// Writting solutions
	extern int write_dual_multipliers(char *workdir, const struct string_buffer *scenarios,
		const struct string_buffer *firststagecols, double **x);
	extern int write_primal_solution(char *workdir, const struct string_buffer *scenarios,
		const struct string_buffer *firststagecols, const double *u, const double *UB);
	
	// Others
	extern double subgradnorm2(const int numscenarios, const double *scenprobabilities,
		const int numfirststagecols, const double *update_scaling, double **v, double *u);
	extern double dual_share(const int k, const int numscenarios, const struct options *opt);
	extern double stepsize(const int k, const double gnorm2, const double LB, const double UB,
		const int numscenarios, const struct options *opt);
	extern int next_job(const int numscenarios, const int numworkers,
		const int numdualworkers, const int numprimalworkers, const double dualshare);
	extern void update_time_stats(int *count, double *average, double *stddev, double new_value);
	extern int time_critical_queue_length(const double av_projection_time, const double stddev_projection_time,
		const int num_arrivals, const double av_arrival_time, const double stddev_arrival_time,
		const int max_primal_workers);
	
	
#endif