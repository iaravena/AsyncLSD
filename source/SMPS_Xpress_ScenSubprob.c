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
 * Functions for manipulating subproblems using Xpress C API
 */

/* SMPS Xpress algorithm header */
#include "SMPS_Xpress.h"

/* Function to convert CORE problem to SCENARIO subproblem */
int transform_CORE_to_SCEN(XPRSprob *problem, struct scen_differences *scen_diff_values)
{
	int i, j, k, l, nrows, ncols;
	XPRSgetintattrib(*problem, XPRS_ROWS, &nrows);
	XPRSgetintattrib(*problem, XPRS_COLS, &ncols);
	
	// Backup CORE coefficient values
	// sorting changes by rows
	if( (*scen_diff_values).numchgcoeff > 0 ){
		int *chgcoeffix;
		chgcoeffix = (int*) malloc((*scen_diff_values).numchgcoeff*sizeof(int));
		for(i = 0; i < (*scen_diff_values).numchgcoeff; i++) *(chgcoeffix + i) = i;
		struct data_frame chgcoeff;
		chgcoeff.numcols = 2;												// data frame with 2 columns
		chgcoeff.numrows = (*scen_diff_values).numchgcoeff;					// I'm not going to acctually copy data, just point to it
		chgcoeff.ptr = (void **) malloc(chgcoeff.numcols*sizeof(void*));
		*(chgcoeff.ptr) = (void *) (*scen_diff_values).chgcoeff_row;		// indexing first by row
		*(chgcoeff.ptr + 1) = (void *) (*scen_diff_values).chgcoeff_col;
		sort_r(chgcoeffix, chgcoeff.numrows, sizeof(int), &compare_rows_idf, (void*) &chgcoeff);
		free(chgcoeff.ptr);
		// retrive CORE values from Xpress
		int nonzeros, retnonzeros, *mstart, *mclind, *mclind_sort;
		struct data_frame mclind_sort_cntxt;
		double *dmatval;
		i = 0;
		while(i < (*scen_diff_values).numchgcoeff){
			// determine range of changes with the same row
			j = 1;
			while( i + j < (*scen_diff_values).numchgcoeff
				&& *((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)) ==
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i + j)) ){
				j++;
			}
			//retrieve CORE values depending on number of changes j
			if(j == 1){
				// get one coefficient at the time
				XPRSgetcoef(*problem,
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)),
					*((*scen_diff_values).chgcoeff_col + *(chgcoeffix + i)),
					((*scen_diff_values).chgcoeff_coreval + *(chgcoeffix + i)));
			} else {
				// get coefficients for an entire row
				XPRSgetrows(*problem, NULL, NULL, NULL, 0, &nonzeros,
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)),
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)));
				mstart = (int *) malloc(2*sizeof(int));					// 1 row requires 2 mstart values
				mclind = (int *) malloc(nonzeros*sizeof(int));
				dmatval = (double *) malloc(nonzeros*sizeof(double));
				if( XPRSgetrows(*problem, mstart, mclind, dmatval, nonzeros, &retnonzeros,
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)),
					*((*scen_diff_values).chgcoeff_row + *(chgcoeffix + i)))
					|| retnonzeros != nonzeros ){
					printf("\nProblem found while retrievieng multiple coefficients from Xpress.\n");
					return(-1);
				}
				// sort column indices for the retrieved coefficients
				mclind_sort = (int *) malloc(nonzeros*sizeof(int));
				for(k = 0; k < nonzeros; k++) *(mclind_sort + k) = k;
				mclind_sort_cntxt.numcols = 1;
				mclind_sort_cntxt.numrows = nonzeros;
				mclind_sort_cntxt.ptr = (void **) malloc(mclind_sort_cntxt.numcols*sizeof(void*));
				*(mclind_sort_cntxt.ptr) = mclind;
				sort_r(mclind_sort, nonzeros, sizeof(int), &compare_rows_idf, (void*) &mclind_sort_cntxt);
				free(mclind_sort_cntxt.ptr);
				// assign relevant values to the restore point
				l = 0;
				for(k = 0; k < j; k++){
					while( *((*scen_diff_values).chgcoeff_col + *(chgcoeffix + i + k)) >
							*(mclind + *(mclind_sort + l))
						&& l < nonzeros - 1) {
						l++;
					}
					if( *((*scen_diff_values).chgcoeff_col + *(chgcoeffix + i + k)) == *(mclind + *(mclind_sort + l)) ){
						*((*scen_diff_values).chgcoeff_coreval + *(chgcoeffix + i + k)) = *(dmatval + *(mclind_sort + l));
					} else {
						*((*scen_diff_values).chgcoeff_coreval + *(chgcoeffix + i + k)) = 0.0;
					}
				}
				free(mstart);
				free(mclind);
				free(dmatval);
				free(mclind_sort);
			}
			i = i + j;
		}
		free(chgcoeffix);
	}
	
	// Backup CORE rhs values
	if( (*scen_diff_values).numchgrhs > 0 ){
		double *CORErhs;
		CORErhs = (double *) malloc(nrows*sizeof(double));
		if( XPRSgetrhs(*problem, CORErhs, 0, nrows-1) != 0 ){
			printf("\nError retrivieng RHS from Xpress.\n");
			return(-1);
		}
		for(i = 0; i < (*scen_diff_values).numchgrhs; i++){
			*((*scen_diff_values).chgrhs_coreval + i) = *(CORErhs + *((*scen_diff_values).chgrhs_row + i));
		}
		free(CORErhs);
	}
	
	// Backup CORE bound values
	if( (*scen_diff_values).numchgbnd > 0 ){
		double *COREubnd, *CORElbnd;
		COREubnd = (double *) malloc(ncols*sizeof(double));
		if( XPRSgetub(*problem, COREubnd, 0, ncols-1) != 0 ){
			printf("\nError retrivieng UB from Xpress.\n");
			return(-1);
		}
		for(i = 0; i < (*scen_diff_values).numchgbnd; i++){
			//if( strncmp((*scen_diff_values).chgbnd_type + i, "U", 1) ){ -> bug! always use strcmp == 0
			if( *((*scen_diff_values).chgbnd_type + i) == 'U' ){
				*((*scen_diff_values).chgbnd_coreval + i) =
					*(COREubnd + *((*scen_diff_values).chgbnd_col + i));
			}
		}
		free(COREubnd);
		CORElbnd = (double *) malloc(ncols*sizeof(double));
		if( XPRSgetlb(*problem, CORElbnd, 0, ncols-1) != 0 ){
			printf("\nError retrivieng LB from Xpress.\n");
			return(-1);
		}
		for(i = 0; i < (*scen_diff_values).numchgbnd; i++){
			//if( strncmp((*scen_diff_values).chgbnd_type + i, "L", 1) ){ -> bug! always use strcmp == 0
			if( *((*scen_diff_values).chgbnd_type + i) == 'L' ){
				*((*scen_diff_values).chgbnd_coreval + i) =
					*(CORElbnd + *((*scen_diff_values).chgbnd_col + i));
			}
		}
		free(CORElbnd);
	}
	
	// Change coefficients to SCEN values
	if( XPRSchgmcoef(*problem, (*scen_diff_values).numchgcoeff, (*scen_diff_values).chgcoeff_row,
		(*scen_diff_values).chgcoeff_col, (*scen_diff_values).chgcoeff_scenval) != 0 ){
		printf("\nError while modifying coefficients in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Change rhs to SCEN values
	if( XPRSchgrhs(*problem, (*scen_diff_values).numchgrhs, (*scen_diff_values).chgrhs_row,
		(*scen_diff_values).chgrhs_scenval) != 0 ){
		printf("\nError while modifying RHS in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Change bounds to SCEN values
	if( XPRSchgbounds(*problem, (*scen_diff_values).numchgbnd, (*scen_diff_values).chgbnd_col,
		(*scen_diff_values).chgbnd_type, (*scen_diff_values).chgbnd_scenval) != 0 ){
		printf("\nError while modifying bounds in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Function to convert SCENARIO subproblem to CORE problem */
int transform_SCEN_to_CORE(XPRSprob *problem, const struct scen_differences *scen_diff_values)
{
	// Change coefficients to CORE values
	if( XPRSchgmcoef(*problem, (*scen_diff_values).numchgcoeff, (*scen_diff_values).chgcoeff_row,
		(*scen_diff_values).chgcoeff_col, (*scen_diff_values).chgcoeff_coreval) != 0 ){
		printf("\nError while modifying coefficients in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Change rhs to CORE values
	if( XPRSchgrhs(*problem, (*scen_diff_values).numchgrhs, (*scen_diff_values).chgrhs_row,
		(*scen_diff_values).chgrhs_coreval) != 0 ){
		printf("\nError while modifying RHS in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Change bounds to CORE values
	if( XPRSchgbounds(*problem, (*scen_diff_values).numchgbnd, (*scen_diff_values).chgbnd_col,
		(*scen_diff_values).chgbnd_type, (*scen_diff_values).chgbnd_coreval) != 0 ){
		printf("\nError while modifying bounds in the matrix held by Xpress.\n");
		return(-1);
	}
	
	// Return success indicator
	return(0);
}

/* Function to add dual multipliers to objective of scenario subproblem */
int load_scenario_multipliers(XPRSprob *scenprob, const int num_firststagecols,
	const int *firststagecols, const double *firststagelb, const double *firststageub,
	const double *scenario_multipliers, double *original_objective)
{
	// Iterators
	int i, j;
	
	// Get current objective values
	int ncolscall = (*(firststagecols + num_firststagecols - 1)) - (*firststagecols) + 1;
	double *obj = (double*) malloc(sizeof(double)*ncolscall);
	XPRSgetobj(*scenprob, obj, *(firststagecols), *(firststagecols + num_firststagecols - 1));
	for(i = 0; i < num_firststagecols; i++){
		j = *(firststagecols + i) - *(firststagecols);
		*(original_objective + i) = *(obj + j);
	}
	free(obj);
	
	// Modify objective on Xpress
	double *new_objective = (double*) malloc(sizeof(double)*num_firststagecols);
	for(i = 0; i < num_firststagecols; i++){
		if( (*(firststageub + i)) - (*(firststagelb + i)) > DBL_EPSILON ){
			*(new_objective + i) = (*(original_objective + i)) +
				1.0/((*(firststageub + i)) - (*(firststagelb + i))) * (*(scenario_multipliers + i));
		} else {
			*(new_objective + i) = (*(original_objective + i)) + (*(scenario_multipliers + i));
		}
	}
	XPRSchgobj(*scenprob, num_firststagecols, firststagecols, new_objective);
	free(new_objective);
	
	// Return success indicator
	return(0);
}

/* Function to restore objective (without multipliers) of scenario subproblem */
int unload_scenario_multipliers(XPRSprob *scenprob, const int num_firststagecols,
	const int *firststagecols, const double *original_objective)
{
	// Modify objective on Xpress
	XPRSchgobj(*scenprob, num_firststagecols, firststagecols, original_objective);
	
	// Return success indicator
	return(0);
}

/* Function to fix first stage variables in scenario subproblem */
int fix_first_stage_vars(XPRSprob *scenprob, const int num_firststagecols,
	const int *firststagecols, double *candidate_sol,
	double *firststagelb, double *firststageub)
{
	// firststagecols should be sorted in ascending order
	// firststagecols, candidate_sol, firststagelb, firststageub should be of lenght num_firststagecols
	
	// Iterators
	int i, j;
	
	// Get current bounds on first stage variables
	int ncolscall = (*(firststagecols + num_firststagecols - 1)) - (*firststagecols) + 1;
	double *dlb = (double*) malloc(sizeof(double)*ncolscall);
	XPRSgetlb(*scenprob, dlb, *(firststagecols), *(firststagecols + num_firststagecols - 1));
	double *dub = (double*) malloc(sizeof(double)*ncolscall);
	XPRSgetub(*scenprob, dub, *(firststagecols), *(firststagecols + num_firststagecols - 1));
	for(i = 0; i < num_firststagecols; i++){
		j = (*(firststagecols + i)) - (*firststagecols);
		*(firststagelb + i) = *(dlb + j);
		*(firststageub + i) = *(dub + j);
	}
	free(dlb);
	free(dub);
	
	// Adapt candidate to current bounds and round integer and binary candidate_sol values
	char *coltype = (char*) malloc(sizeof(char)*ncolscall);
	XPRSgetcoltype(*scenprob, coltype, *(firststagecols), *(firststagecols + num_firststagecols - 1));
	for(i = 0; i < num_firststagecols; i++){
		j = (*(firststagecols + i)) - (*firststagecols);
		if( *(coltype + j) == 'B' || *(coltype + j) == 'I' ){
			*(candidate_sol + i) = rint(*(candidate_sol + i));
		}
	}
	free(coltype);
	
	// Fix first stage variables 
	char *bndtype = (char*) malloc(sizeof(char)*num_firststagecols);
	for(i = 0; i < num_firststagecols; i++) *(bndtype + i) = 'B';
	XPRSchgbounds(*scenprob, num_firststagecols, firststagecols, bndtype, candidate_sol);
	free(bndtype);
	
	// Return success indicator
	return(0);
}

/* Function to unfix first stage variables in scenario subproblem */
int unfix_first_stage_vars(XPRSprob *scenprob, const int num_firststagecols, const int *firststagecols,
	const double *firststagelb, const double *firststageub)
{
	// firststagecols should be sorted in ascending order
	// firststagecols, firststagelb, firststageub should be of lenght num_firststagecols
	
	// Iterators and bound type variable
	int i;
	char *bndtype;
	
	// Restore lower bounds
	bndtype = (char*) malloc(sizeof(char)*num_firststagecols);
	for(i = 0; i < num_firststagecols; i++) *(bndtype + i) = 'L';
	XPRSchgbounds(*scenprob, num_firststagecols, firststagecols, bndtype, firststagelb);
	free(bndtype);
	
	// Restore upper bounds
	bndtype = (char*) malloc(sizeof(char)*num_firststagecols);
	for(i = 0; i < num_firststagecols; i++) *(bndtype + i) = 'U';
	XPRSchgbounds(*scenprob, num_firststagecols, firststagecols, bndtype, firststageub);
	free(bndtype);
	
	// Return success indicator
	return(0);
}

/* Evaluate sequential period relaxation of scenario subproblem */
int solve_sequential_period_relaxation(const XPRSprob *scenproblem,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const int num_delayed_rows, const int delayed_rows_start,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const int firststageind, const int num_firststgcols, const int *firststgcols,
	double *obj_lb, double *firststgsol, const struct options *opt)
{
	// Iterators
	int i, t;
	
	// Initialize objective lower bound and solution
	*obj_lb = 0;
	if( (*opt).period_aggregation_policy == PERIOD_POLICY_MAX ){
		double *lb = (double*) malloc(sizeof(double) * (*colnames).num_elems);
		XPRSgetlb(*scenproblem, lb, *(firststgcols), *(firststgcols + num_firststgcols - 1));
		for(i = 0; i < num_firststgcols; i++)
			*(firststgsol + i) = *(lb + (*(firststgcols + i)) - *firststgcols);
		free(lb);
	}else if( (*opt).period_aggregation_policy == PERIOD_POLICY_MEAN ){
		for(i = 0; i < num_firststgcols; i++) *(firststgsol + i) = 0;
	}else if( (*opt).period_aggregation_policy == PERIOD_POLICY_MIN ){
		double *ub = (double*) malloc(sizeof(double) * (*colnames).num_elems);
		XPRSgetub(*scenproblem, ub, *(firststgcols), *(firststgcols + num_firststgcols - 1));
		for(i = 0; i < num_firststgcols; i++)
			*(firststgsol + i) = *(ub + (*(firststgcols + i)) - *firststgcols);
		free(ub);
	}else{
		printf("\nUnrecognized period aggregation policy.\n");
		return(-1);
	}
	
	// Determine multiplicity of first stage variables in subsequent periods
	int *firststagecols_mult = malloc(sizeof(int)*num_firststgcols);
	get_firststagecols_multiplicity(scenproblem, rownames, colnames, timestages,
		rowtstages, coltstages, firststageind, num_firststgcols, firststgcols, firststagecols_mult);
	
	// Compute contribution of variables appearing exclusively in the first stage
	double objcoeff = 0, colbnd;
	for(i = 0; i < num_firststgcols; i++){
		if( *(firststagecols_mult + i) == 0
			&& XPRSgetobj(*scenproblem, &objcoeff, *(firststagecols_mult + i), *(firststagecols_mult + i)) == 0
			&& objcoeff != 0){	// NOTE: I am minimizing
			if(objcoeff > 0){
				XPRSgetlb(*scenproblem, &colbnd, *(firststagecols_mult + i), *(firststagecols_mult + i));
			} else {
				XPRSgetub(*scenproblem, &colbnd, *(firststagecols_mult + i), *(firststagecols_mult + i));;
			}
			*obj_lb = *obj_lb + objcoeff*colbnd;
			*(firststgsol + i) = colbnd;
		}
	}
	
	// Period problem attributes (debbuging only)
	//int ncols, nrows, nelems, nglobal;
	//printf("\n\nPeriod subproblems:\n"
	//	"Ts\tCols\tRows\tElems\tGlob\tOpt\n");
	
	// Loop over all periods but first stage 
	XPRSprob scenperiodprob;
	int numincfscols, *incfscols, *incfscols_ind, nscenpercols, mipstatus;
	double *scenperiodsol, scenperobj;
	struct timespec startt, endt;
	incfscols = (int*) malloc(sizeof(int)*num_firststgcols);
	incfscols_ind = (int*) malloc(sizeof(int)*num_firststgcols);
	for(t = 0; t < (*timestages).num_elems; t++){
		if( t != firststageind ){											// omitting T0, first stage
			
			// Create period subproblem object in Xpress
			XPRScreateprob(&scenperiodprob);
			#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
				XPRSaddcbmessage(scenperiodprob, Message, NULL, 0);			// defines callback function for printing log from Xpress (usefull only in Windows)
			#endif
			
			// load xpress parameters
			{
				int j;
				for(j = 0; j < (*opt).periodprob_options.num_intoptions; j++){
					XPRSsetintcontrol(scenperiodprob, *((*opt).periodprob_options.intoptions + j),
						*((*opt).periodprob_options.intoptions_value + j));
				}
				for(j = 0; j < (*opt).periodprob_options.num_dbloptions; j++){
					XPRSsetdblcontrol(scenperiodprob, *((*opt).periodprob_options.dbloptions + j),
						*((*opt).periodprob_options.dbloptions_value + j));
				}
			}
			
			// override verbosity level
			if((*opt).verbosity > VERB_WORKER_FUNCS){
				XPRSsetintcontrol(scenperiodprob, XPRS_OUTPUTLOG, 1);
			} else {
				XPRSsetintcontrol(scenperiodprob, XPRS_OUTPUTLOG, 4);
			}
			
			// Load period subproblem parameters
			current_utc_time(&startt);
			if((*opt).verbosity >= VERB_WORKER_FUNCS)
				printf("\nLoading subproblem for period %d ...\n", t);
			if( load_params_period_subprob(scenproblem, "period_problem", rownames, colnames,
				num_delayed_rows, delayed_rows_start,
				timestages, rowtstages, coltstages, firststageind, t,
				num_firststgcols, firststgcols, firststagecols_mult,
				&scenperiodprob, &numincfscols, incfscols, incfscols_ind) != 0){
				printf("\nError while loading subproblem for period %d.\n", t);
				return(-1);
			}
			current_utc_time(&endt);
			if((*opt).verbosity >= VERB_WORKER_FUNCS)
				printf("\nSubproblem for period %d loaded in %.2f seconds.\n", t,
					((double) (endt.tv_sec - startt.tv_sec)) +
					((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC);
			
			// Get problem parameters (debugging only)
			//XPRSgetintattrib(scenperiodprob, XPRS_ORIGINALCOLS, &ncols);
			//XPRSgetintattrib(scenperiodprob, XPRS_ORIGINALROWS, &nrows);
			//XPRSgetintattrib(scenperiodprob, XPRS_ELEMS, &nelems);
			//XPRSgetintattrib(scenperiodprob, XPRS_ORIGINALMIPENTS, &nglobal);
			
			// Solve period subproblem
			current_utc_time(&startt);
			if((*opt).verbosity >= VERB_WORKER_FUNCS)
				printf("\nSolving subproblem for period %d ...\n", t);
			XPRSmipoptimize(scenperiodprob, (*opt).periodprob_options.lp_method);
			XPRSgetintattrib(scenperiodprob, XPRS_MIPSTATUS, &mipstatus);
			if( !(mipstatus == XPRS_MIP_OPTIMAL || mipstatus == XPRS_MIP_SOLUTION) ){
				printf("\nXpress failed at solving period relaxation for period %s."
					"\nProblem status (MIPSTATUS):\t", (*timestages).buffer + ((*timestages).atom_len+1)*t);
				switch(mipstatus){
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
				return(-1);
			}
			current_utc_time(&endt);
			if((*opt).verbosity >= VERB_WORKER_FUNCS)
				printf("\nSubproblem for period %d solved in %.2f seconds.\n", t,
					((double) (endt.tv_sec - startt.tv_sec)) +
					((double) (endt.tv_nsec - startt.tv_nsec))*NANOSEC2SEC);
			
			// Store bound contribution and first stage solution
			XPRSgetintattrib(scenperiodprob, XPRS_COLS, &nscenpercols);
			scenperiodsol = (double*) malloc(sizeof(double)*nscenpercols);
			XPRSgetmipsol(scenperiodprob, scenperiodsol, NULL);
			for(i = 0; i < numincfscols; i++){
				if( (*opt).period_aggregation_policy == PERIOD_POLICY_MAX ){
					*(firststgsol + *(incfscols_ind + i)) =
						fmax((*(firststgsol + *(incfscols_ind + i))), *(scenperiodsol + *(incfscols + i)));
				}else if( (*opt).period_aggregation_policy == PERIOD_POLICY_MEAN ){
					*(firststgsol + *(incfscols_ind + i)) +=
						*(scenperiodsol + *(incfscols + i))/((double) *(firststagecols_mult + *(incfscols_ind + i)));
				}else{ // PERIOD_POLICY_MEAN
					*(firststgsol + *(incfscols_ind + i)) =
						fmin((*(firststgsol + *(incfscols_ind + i))), *(scenperiodsol + *(incfscols + i)));
				}
			}
			free(scenperiodsol);
			XPRSgetdblattrib(scenperiodprob, XPRS_BESTBOUND, &scenperobj);
			*obj_lb = (*obj_lb) + scenperobj;
			
			// Delete period subproblem from Xpress
			XPRSdestroyprob(scenperiodprob);
			
			// Write problem attributes (debugging only)
			//printf("%d\t%d\t%d\t%d\t%d\t%f\n", t, ncols, nrows, nelems, nglobal, scenperobj);
		}
	}
	
	// Get problem parameters (debugging only)
	//XPRSgetintattrib(*scenproblem, XPRS_ORIGINALCOLS, &ncols);
	//XPRSgetintattrib(*scenproblem, XPRS_ORIGINALROWS, &nrows);
	//XPRSgetintattrib(*scenproblem, XPRS_ELEMS, &nelems);
	//XPRSgetintattrib(*scenproblem, XPRS_ORIGINALMIPENTS, &nglobal);
	//printf("\n\nScenario problem:\n"
	//	"Cols\tRows\tElems\tGlob\n"
	//	"%d\t%d\t%d\t%d\n\n", ncols, nrows, nelems, nglobal);
	
	// Free heap
	free(firststagecols_mult);
	free(incfscols);
	free(incfscols_ind);
	
	// Return success indicator
	return(0);
}
