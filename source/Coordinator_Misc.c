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

/* Write dual multpliers to CSV */
int write_dual_multipliers(char *workdir, const struct string_buffer *scenarios,
	const struct string_buffer *firststagecols, double **x)
{
	char multiplierfile[200];
	FILE *out;
	int i, j;
	
	strcpy(multiplierfile, workdir);
	strcat(multiplierfile, SYSTEM_SLASH);
	strcat(multiplierfile, "dual_multipliers.csv");
	out = fopen(multiplierfile, "w");
	fprintf(out, "j");
	for(i = 0; i < (*scenarios).num_elems; i++)
		fprintf(out, ",%s", (*scenarios).buffer + ((*scenarios).atom_len + 1) * i);
	fprintf(out, "\n");
	for(j = 0; j < (*firststagecols).num_elems; j++){
		fprintf(out, "%s", (*firststagecols).buffer + ((*firststagecols).atom_len + 1) * j);
		for(i = 0; i < (*scenarios).num_elems; i++){
			fprintf(out, ",%f", *((*(x + i)) + j));
		}
		fprintf(out, "\n");
	}
	fclose(out);
	
	return(0);
}

/* Write incumbent to CSV files */
int write_primal_solution(char *workdir, const struct string_buffer *scenarios,
	const struct string_buffer *firststagecols, const double *u, const double *UB)
{
	char incumbentfile[200], incumbentboundsfile[200];
	FILE *out;
	int i;
	
	// incumbent file
	strcpy(incumbentfile, workdir);
	strcat(incumbentfile, SYSTEM_SLASH);
	strcat(incumbentfile, "primal_incumbent.csv");
	out = fopen(incumbentfile, "w");
	fprintf(out, "Column,Value\n");
	for(i = 0; i < (*firststagecols).num_elems; i++){
		fprintf(out, "%s,%f\n",
			(*firststagecols).buffer + ((*firststagecols).atom_len + 1) * i,
			*(u + i));
	}
	fclose(out);
	
	// incumbent bounds file
	strcpy(incumbentboundsfile, workdir);
	strcat(incumbentboundsfile, SYSTEM_SLASH);
	strcat(incumbentboundsfile, "primal_incumbent_UB.csv");
	out = fopen(incumbentboundsfile, "w");
	fprintf(out, "Scenario,UB\n");
	for(i = 0; i < (*scenarios).num_elems; i++){
		fprintf(out, "%s,%f\n",
		(*scenarios).buffer + ((*scenarios).atom_len + 1) * i,
		*(UB + i));
	}
	fclose(out);
	
	return(0);
}

/* Compute subgradient norm */
double subgradnorm2(const int numscenarios, const double *scenprobabilities,
	const int numfirststagecols, const double *update_scaling, double **v, double *u)
{
	int i, j;
	double gnorm2 = 0;
	
	for(i = 0; i < numscenarios; i++){
		for(j = 0; j < numfirststagecols; j++){
			if( *(update_scaling + j) > NO_UPDATE_SCALING ){
				gnorm2 = gnorm2 + pow((*(scenprobabilities + i)) * (*(update_scaling + j)) *
					((*((*(v + i)) + j)) - *(u + j)), 2.0);
			}else{
				gnorm2 = gnorm2 + pow((*(scenprobabilities + i)) * ((*((*(v + i)) + j)) - *(u + j)), 2.0);
			}
			//gnorm2 += pow((*((*(v + i)) + j)) - (*(u + j)), 2);
		}
	}
	
	return(gnorm2);
}

/* Dual share dynamic computation */
double dual_share(const int k, const int numscenarios, const struct options *opt)
{
	int max_CD_iterations = numscenarios * (*opt).max_iterations - 1;
	double dsmin = tanh(-5.0 * (*opt).dual_share_midpoint);
	double dsmax = tanh(5.0 * (1 - (*opt).dual_share_midpoint));
	double ds = ((*opt).dual_share_end - (*opt).dual_share_start)/(dsmax - dsmin) *
		(tanh(5.0/max_CD_iterations*k - 5.0*(*opt).dual_share_midpoint)-dsmin) + (*opt).dual_share_start;
	return(ds);
}

/* Compute stepsize */
double stepsize(const int k, const double gnorm2, const double LB, const double UB,
	const int numscenarios, const struct options *opt)
{
	double lambda=0;
	if( (*opt).dual_stepsize == DUAL_STEPSIZE_CONSTANT ){
		lambda = (*opt).dual_stepsize_p;
	}else if( (*opt).dual_stepsize == DUAL_STEPSIZE_DIMINISHING ){
		lambda = (*opt).dual_stepsize_p/pow(1 + (*opt).dual_stepsize_r/numscenarios*k, (*opt).dual_stepsize_q);
	}else if( (*opt).dual_stepsize == DUAL_STEPSIZE_POLYAK ){
		lambda = (*opt).dual_stepsize_p/pow(1 + (*opt).dual_stepsize_r/numscenarios*k, (*opt).dual_stepsize_q)*
			fmin( (*opt).dual_stepsize_xi*fabs(LB), UB-LB )/fmax((*opt).dual_stepsize_sigma, gnorm2);
	}
	return(lambda);
}

/* Routine to determine next job to run */
int next_job(const int numscenarios, const int numworkers,
	const int numdualworkers, const int numprimalworkers, const double dualshare)
{
	int target_dualworkers = (int) fmin((double) numscenarios, floor(((double) numworkers)*dualshare));
	int target_primalworkers = numworkers - 1 - target_dualworkers;
	if(numdualworkers < target_dualworkers){
		return(NEXT_COMES_DUAL);
	}else if(numprimalworkers < target_primalworkers){
		return(NEXT_COMES_PRIMAL);
	}else{
		double dualremainder = dualshare - ((double) target_dualworkers)/((double) numworkers);
		if( rbernoulli(dualremainder) ){
			return(NEXT_COMES_DUAL);
		}else{
			return(NEXT_COMES_PRIMAL);
		}
	}
}

/* Update time statistics */
void update_time_stats(int *count, double *average, double *stddev, double new_value)
{
	// compute sum(x) and sum(x^2)
	double sum, sum2;
	sum = (*average) * (*count) + new_value;
	sum2 = pow((*stddev), 2.0) * ((*count) - 1) + (*count) * pow(*average, 2.0) + pow(new_value, 2.0);
	
	// recompute average and standard deviation
	*count += 1;
	*average = sum/((double) *count);
	*stddev = pow((sum2 - (*count) * pow(*average, 2.0))/(((double) *count) - 1.0), 0.5);
}

/* Time critical primal queue length */
int time_critical_queue_length(const double av_projection_time, const double stddev_projection_time,
	const int num_arrivals, const double av_arrival_time, const double stddev_arrival_time,
	const int max_primal_workers)
{
	double projection_time_UB95 = av_projection_time + QNORM95*stddev_projection_time;
	double arrival_time_LB95 = fmax(0.0, av_arrival_time -
		QNORM95/pow((double) num_arrivals, 0.5)*stddev_arrival_time);
	int queue_length = (int) fmin((double) MAX_LENGTH_CRITICAL_QUEUE*max_primal_workers,
		ceil(projection_time_UB95/arrival_time_LB95));
	//printf("\nCritical queue length:\t%f\t%f\t%d\n", projection_time_UB95,
	//	arrival_time_LB95, queue_length);
	return(queue_length);
}