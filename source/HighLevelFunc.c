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
 * High level functions
 * IMPORTANT: these functions include calls to malloc that are not inmediately released, do not forget to clean!
 */

/* Asynchronous algorithm header */
#include "AsyncHeader.h"

/* Wrapper for reading CORE file, TIME file and general information on scenarios
 * (passing pointers to pointers only to be able to malloc and free) */
int read_SMPS(const char *workdir, const char *SMPSrootname,
	XPRSprob *problem, struct string_buffer *rownames, struct string_buffer *colnames,
	struct string_buffer *timestages, int **rowtstages, int **coltstages,
	struct string_buffer *scenarios, struct string_buffer *scenparent, int **scenbranchtstage,
	int **scenstartline, int **scenendline, double **scenprobability,
	int *firststageind, int *numfirststagecols, int **firststagecols,
	double **fstg_lbnds, double **fstg_ubnds, const struct options *opt)
{
	// Set output from Xpress to the verbosity level
	int xpress_output_log;
	XPRSgetintcontrol(*problem, XPRS_OUTPUTLOG, &xpress_output_log);
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 4);
	}
	
	// Load CORE problem to Xpress
	char mpsname[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nLoading CORE file to Xpress ...\n");
	strcpy(mpsname, workdir);
	strcat(strcat(strcat(mpsname, SYSTEM_SLASH), SMPSrootname), ".mps");
	XPRSreadprob(*problem, mpsname, "");
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("Problem loaded in Xpress.\n");
	
	// Read problem row and column names
	int nl;
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nGathering row and column names from Xpress ...");
	XPRSgetintattrib(*problem, XPRS_ROWS, &((*rownames).num_elems));
	XPRSgetintattrib(*problem, XPRS_COLS, &((*colnames).num_elems));
	XPRSgetintattrib(*problem, XPRS_NAMELENGTH, &nl);
	(*rownames).atom_len = uXPRS_NAMELENGTH_DEF*nl;
	(*colnames).atom_len = uXPRS_NAMELENGTH_DEF*nl;
	(*rownames).buffer = (char *) malloc(sizeof(char)*((*rownames).atom_len+1)*(*rownames).num_elems);
	(*colnames).buffer = (char *) malloc(sizeof(char)*((*colnames).atom_len+1)*(*colnames).num_elems);
	XPRSgetnames(*problem,1,(*rownames).buffer,0,(*rownames).num_elems-1);
	XPRSgetnames(*problem,2,(*colnames).buffer,0,(*colnames).num_elems-1);
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nRows and column names available to main program.\n");
	
	// Read TIME file
	char timename[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nParsing TIME file ...");
	strcpy(timename, workdir);
	strcat(strcat(strcat(timename, SYSTEM_SLASH), SMPSrootname), ".tim");
	(*timestages).atom_len = SMPS_MAXFIELDLEN;
	(*timestages).num_elems = 0;
	(*timestages).buffer = (char *) malloc(sizeof(char)*((*timestages).atom_len+1)*SMPS_MAXSTAGES);
	*rowtstages = (int *) calloc((*rownames).num_elems, sizeof(int));
	*coltstages = (int *) calloc((*colnames).num_elems, sizeof(int));
	if( read_TIME_file(SMPSrootname, timename, rownames, colnames, timestages, *rowtstages, *coltstages) ){
		printf("\nError while reading TIME file.\n");
		return(-1);
		}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\n%d time stages detected. All rows and cols correctly assigned to time stages.\n",
			(*timestages).num_elems);
	
	// Read general scenario info from STOCH file
	char stochname[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nPreprocessing STOCH file:\n(i) Getting number and name length of scenarios ...");
	strcpy(stochname, workdir);
	strcat(strcat(strcat(stochname, SYSTEM_SLASH), SMPSrootname), ".sto");
	
	// get scenario set size for creating buffers
	if( get_scenario_sizes(SMPSrootname, stochname, &((*scenarios).atom_len), &((*scenarios).num_elems)) ){
		printf("\nError while reading the character length and number of scenarios.\n");
		return(-1);
		}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\n%d scenarios with names up to %d characters long detected in STOCH file.\n",
			(*scenarios).num_elems, (*scenarios).atom_len);
	
	// read scenario metadata into buffers
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("(ii) Allocating space and reading scenario metadata ...");
	(*scenarios).buffer = (char *) malloc(sizeof(char)*((*scenarios).atom_len+1)*(*scenarios).num_elems);
	(*scenparent).atom_len = (int) fmax((*scenarios).atom_len, strlen(SMPS_ROOTNODENAME));
	(*scenparent).num_elems = (*scenarios).num_elems;
	(*scenparent).buffer = (char *) malloc(sizeof(char)*((*scenparent).atom_len+1)*(*scenparent).num_elems);
	*scenbranchtstage = (int *) calloc((*scenarios).num_elems, sizeof(int));
	*scenstartline = (int *) calloc((*scenarios).num_elems, sizeof(int));
	*scenendline = (int *) calloc((*scenarios).num_elems, sizeof(int));
	*scenprobability = (double *) calloc((*scenarios).num_elems, sizeof(double));
	if( get_scenario_metadata(SMPSrootname, stochname, timestages,
		scenarios, scenparent, *scenbranchtstage, *scenprobability, *scenstartline, *scenendline) ){
		printf("\nError while reading scenarios' metadata.\n");
		return(-1);
		}
	
	// verify that the probabilities add up to 1
	{
		int i;
		double totalprob = 0;
		for(i = 0; i < (*scenarios).num_elems; i++){
			totalprob += *((*scenprobability) + i);
		}
		if( totalprob < 1.0){
			//message only shown by coordinator
			//printf("\nWARNING: Scenario probalities are %e less than 1. They will be rescaled up to 1.\n", 1.0-totalprob);
			for(i = 0; i < (*scenarios).num_elems; i++){
				*((*scenprobability) + i) /= totalprob;
			}
		}
	}
	
	// print scenario data summary
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		int i, minchg, maxchg;
		float minprob, maxprob;
		minchg = *((*scenendline) + (*scenarios).num_elems - 1);
		maxchg = 0;
		minprob = 1;
		maxprob = 0;
		for(i = 0; i < (*scenarios).num_elems; i++){
			minchg = (int) fmin(minchg, (*((*scenendline) + i)) - (*((*scenstartline) + i)));
			maxchg = (int) fmax(maxchg, (*((*scenendline) + i)) - (*((*scenstartline) + i)));
			minprob = (float) fmin(minprob, *((*scenprobability) + i));
			maxprob = (float) fmax(maxprob, *((*scenprobability) + i));
		}
		printf("\nScenario metadata for %d scenarios read from STOCH file."
			"\nNumber of changes to CORE (min - max):\t%d - %d"
			"\nProbabilities (min - max):\t\t%.2f%% - %.2f%%\n",
			(*scenarios).num_elems, minchg, maxchg, 100*minprob, 100*maxprob);
	}
	
	// Get first stage index and determine which columns are first stage columns
	char firststagefilled[(*timestages).atom_len+1];
	add_spaces(SMPS_FIRSTSTAGENAME, (*timestages).atom_len+1, (*timestages).atom_len, firststagefilled);
	*firststageind = get_string_index(timestages, &firststagefilled[0]);
	if(*firststageind < 0){
		printf("\nFirst stage '" SMPS_FIRSTSTAGENAME "' not found among time stages.\n");
		return(-1);
	}
	*firststagecols = (int*) malloc(sizeof(int)*(*colnames).num_elems);
	select_cols_and_rows(problem, rownames, colnames, 0, 0, timestages, *rowtstages, *coltstages,
		1, firststageind, 0, NULL, numfirststagecols, *firststagecols,
		NULL, NULL, NULL);
	*firststagecols = (int*) realloc(*firststagecols, sizeof(int) * *numfirststagecols);
	
	// Get first stage cols bounds from CORE problem
	*fstg_lbnds = (double*) malloc(sizeof(double) * *numfirststagecols);
	*fstg_ubnds = (double*) malloc(sizeof(double) * *numfirststagecols);
	{
		int i;
		double *bndall = (double*) malloc(sizeof(double)*(*colnames).num_elems);
		XPRSgetlb(*problem, bndall, 0, (*colnames).num_elems-1);
		for(i = 0; i < *numfirststagecols; i++){
			*((*fstg_lbnds) + i) = *(bndall + *((*firststagecols) + i));
		}
		XPRSgetub(*problem, bndall, 0, (*colnames).num_elems-1);
		for(i = 0; i < *numfirststagecols; i++){
			*((*fstg_ubnds) + i) = *(bndall + *((*firststagecols) + i));
		}
		free(bndall);
	}
	
	// Get first stage colnames
	struct string_buffer firststagecolnames;
	firststagecolnames.num_elems = *numfirststagecols;
	firststagecolnames.atom_len = (*colnames).atom_len;
	firststagecolnames.buffer = (char*) malloc(sizeof(char)*
		(firststagecolnames.atom_len+1)*firststagecolnames.num_elems);
	{
		int i;
		for(i = 0; i < firststagecolnames.num_elems; i++){
			insert_in_buffer(&firststagecolnames, i,
				(*colnames).buffer+((*colnames).atom_len+1) * *((*firststagecols) + i));
		}
	}
	
	// Gather scenario bounds from STOCH file, scenario by scenario
	{
		int i;
		for(i = 0; i < (*scenarios).num_elems; i++){
			if( get_scenario_col_bounds(SMPSrootname, stochname,
				(*scenarios).buffer + ((*scenarios).atom_len + 1) * i,
				scenarios, scenparent, *scenstartline, *scenendline,
				&firststagecolnames, *fstg_lbnds, *fstg_ubnds) != 0 ){
				printf("\nProblem reading outer bounds for first stage variables at scenario %s.",
					(*scenarios).buffer + ((*scenarios).atom_len + 1) * i);
				return(-1);
			}
		}
	}
	
	// Free first stage column names
	free(firststagecolnames.buffer);
	
	// Check that all first stage variables are bounded
	{
		int i;
		for(i = 0; i < firststagecolnames.num_elems; i++){
			if( *((*fstg_lbnds) + i) <= XPRS_MINUSINFINITY || *((*fstg_ubnds) + i) >= XPRS_PLUSINFINITY ){
				printf("\nUnbounded column %s detected in first stage matrix.\n",
					(*colnames).buffer + ((*colnames).atom_len + 1) * *((*firststagecols) + i) );
				return(-1);
			}
		}
	}
	
	// Put cross-scenario consistent bounds in the CORE problem
	{
		int i;
		char *bndtype = (char*) malloc(sizeof(char) * *numfirststagecols);
		for(i = 0; i < *numfirststagecols; i++) *(bndtype + i) = 'L';
		XPRSchgbounds(*problem, *numfirststagecols, *firststagecols, bndtype, *fstg_lbnds);
		for(i = 0; i < *numfirststagecols; i++) *(bndtype + i) = 'U';
		XPRSchgbounds(*problem, *numfirststagecols, *firststagecols, bndtype, *fstg_ubnds);
		free(bndtype);
	}
	
	// Restore output from Xpress to the value outside this function
	XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, xpress_output_log);
	
	// Return success indicator
	return(0);
}

/* Wrapper for cleaning malloc calls inside read_SMPS
 * (passing pointers to pointers only to be able to malloc and free) */
int free_heap_from_read_SMPS(
	XPRSprob *problem, struct string_buffer *rownames, struct string_buffer *colnames,
	struct string_buffer *timestages, int **rowtstages, int **coltstages,
	struct string_buffer *scenarios, struct string_buffer *scenparent, int **scenbranchtstage,
	int **scenstartline, int **scenendline, double **scenprobability,
	int **firststagecols, double **fstg_lbnds, double **fstg_ubnds, const struct options *opt)
{
	// Free memory allocated on the HEAP
	free((*rownames).buffer);
	free((*colnames).buffer);
	free((*timestages).buffer);
	free(*rowtstages);
	free(*coltstages);
	free((*scenarios).buffer);
	free((*scenparent).buffer);
	free(*scenbranchtstage);
	free(*scenstartline);
	free(*scenendline);
	free(*scenprobability);
	free(*firststagecols);
	free(*fstg_lbnds);
	free(*fstg_ubnds);
	
	// Return success indicator
	return(0);
}

/* Wrapper for loading a scenario subproblem */
int load_scenario_subproblem(const char *workdir, const char *SMPSrootname, const int scenindex,
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *problem, struct scen_differences *change_and_restoration_point, const struct options *opt)
{
	// Set output from Xpress to the verbosity level
	int xpress_output_log;
	XPRSgetintcontrol(*problem, XPRS_OUTPUTLOG, &xpress_output_log);
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 4);
	}
	
	// Read differences between SCEN and CORE
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nReading STOCH file for scenario %s ...",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*scenindex);
	
	char stochname[200];
	strcpy(stochname, workdir);
	strcat(strcat(strcat(stochname, SYSTEM_SLASH), SMPSrootname), ".sto");
	
	// bounds on the number of changes
	(*change_and_restoration_point).numchgcoeff = (int) fmin((*rownames).num_elems*(*colnames).num_elems,
		(*(scenendline + scenindex)) - (*(scenstartline + scenindex)));
	(*change_and_restoration_point).numchgrhs = (int) fmin((*rownames).num_elems,
		(*(scenendline + scenindex)) - (*(scenstartline + scenindex)));
	(*change_and_restoration_point).numchgbnd = (int) fmin((*colnames).num_elems,
		(*(scenendline + scenindex)) - (*(scenstartline + scenindex)));
	
	// allocate buffers for recording scenario values
	(*change_and_restoration_point).chgcoeff_col =
		(int *) calloc((*change_and_restoration_point).numchgcoeff, sizeof(int));
	(*change_and_restoration_point).chgcoeff_row =
		(int *) calloc((*change_and_restoration_point).numchgcoeff, sizeof(int));
	(*change_and_restoration_point).chgcoeff_scenval =
		(double *) calloc((*change_and_restoration_point).numchgcoeff, sizeof(double));
	(*change_and_restoration_point).chgrhs_row =
		(int *) calloc((*change_and_restoration_point).numchgrhs, sizeof(int));
	(*change_and_restoration_point).chgrhs_scenval =
		(double *) calloc((*change_and_restoration_point).numchgrhs, sizeof(double));
	(*change_and_restoration_point).chgbnd_col =
		(int *) calloc((*change_and_restoration_point).numchgbnd, sizeof(int));
	(*change_and_restoration_point).chgbnd_type =
		(char *) calloc((*change_and_restoration_point).numchgbnd, sizeof(char));
	(*change_and_restoration_point).chgbnd_scenval =
		(double *) calloc((*change_and_restoration_point).numchgbnd, sizeof(double));
	
	// read scenario differences
	if( *(scenendline + scenindex) - *(scenstartline + scenindex) > 0 ){
		if( read_scen_differences(SMPSrootname, &stochname[0],
				(*scenarios).buffer + ((*scenarios).atom_len+1)*scenindex, scenarios, scenparent,
				scenstartline, rownames, colnames, change_and_restoration_point) ){
			printf("\nError while reading data for scenario for scenario %s.\n",
				(*scenarios).buffer + ((*scenarios).atom_len+1)*scenindex);
			return(-1);
		}
	} else {
		(*change_and_restoration_point).numchgcoeff = 0;
		(*change_and_restoration_point).numchgrhs = 0;
		(*change_and_restoration_point).numchgbnd = 0;
	}
	
	// shrink buffers for SCEN values to their current size
	(*change_and_restoration_point).chgcoeff_col =
		(int *) realloc((void *) (*change_and_restoration_point).chgcoeff_col,
			sizeof(int)*(*change_and_restoration_point).numchgcoeff);
	(*change_and_restoration_point).chgcoeff_row =
		(int *) realloc((void *) (*change_and_restoration_point).chgcoeff_row,
			sizeof(int)*(*change_and_restoration_point).numchgcoeff);
	(*change_and_restoration_point).chgcoeff_scenval =
		(double *) realloc((void *) (*change_and_restoration_point).chgcoeff_scenval,
			sizeof(double)*(*change_and_restoration_point).numchgcoeff);
	(*change_and_restoration_point).chgrhs_row =
		(int *) realloc((void *) (*change_and_restoration_point).chgrhs_row,
			sizeof(int)*(*change_and_restoration_point).numchgrhs);
	(*change_and_restoration_point).chgrhs_scenval =
		(double *) realloc((void *) (*change_and_restoration_point).chgrhs_scenval,
			sizeof(double)*(*change_and_restoration_point).numchgrhs);
	(*change_and_restoration_point).chgbnd_col =
		(int *) realloc((void *) (*change_and_restoration_point).chgbnd_col,
			sizeof(int)*(*change_and_restoration_point).numchgbnd);
	(*change_and_restoration_point).chgbnd_type =
		(char *) realloc((void *) (*change_and_restoration_point).chgbnd_type,
			sizeof(char)*(*change_and_restoration_point).numchgbnd);
	(*change_and_restoration_point).chgbnd_scenval =
		(double *) realloc((void *) (*change_and_restoration_point).chgbnd_scenval,
			sizeof(double)*(*change_and_restoration_point).numchgbnd);
	
	// intermediate success message
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nRead differences for scenario %s."
			"\n# COEFF differences:\t%d"
			"\n# RHS differences:\t%d"
			"\n# BND differences:\t%d\n",
			(*scenarios).buffer + ((*scenarios).atom_len+1)*scenindex,
			(*change_and_restoration_point).numchgcoeff,
			(*change_and_restoration_point).numchgrhs,
			(*change_and_restoration_point).numchgbnd);
	
	// Transform CORE problem into scenario subproblem
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nTransforming CORE problem into scenario subproblem ...");
	
	// allocate buffers for CORE values
	(*change_and_restoration_point).chgcoeff_coreval =
		(double *) calloc((*change_and_restoration_point).numchgcoeff, sizeof(double));
	(*change_and_restoration_point).chgrhs_coreval =
		(double *) calloc((*change_and_restoration_point).numchgrhs, sizeof(double));
	(*change_and_restoration_point).chgbnd_coreval =
		(double *) calloc((*change_and_restoration_point).numchgbnd, sizeof(double));
	
	// perform transformation and backup
	if( transform_CORE_to_SCEN(problem, change_and_restoration_point) != 0){
		printf("\nError while generating scenario subporblem.\n");
		return(-1);
	}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nScenario subproblem ready.\n");
	
	// Restore output from Xpress to the value outside this function
	XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, xpress_output_log);
	
	// Return success indicator
	return(0);
}

/* Wrapper for restoring the CORE problem */
int restore_core_problem(
	const struct string_buffer *rownames, const struct string_buffer *colnames,
	const struct string_buffer *timestages, const int *rowtstages, const int *coltstages,
	const struct string_buffer *scenarios,  const struct string_buffer *scenparent, const int *scenbranchtstage,
	const int *scenstartline, const int *scenendline, const double *scenprobability,
	XPRSprob *problem, struct scen_differences *change_and_restoration_point, const struct options *opt)
{
	// Set output from Xpress to the verbosity level
	int xpress_output_log;
	XPRSgetintcontrol(*problem, XPRS_OUTPUTLOG, &xpress_output_log);
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, 4);
	}
	
	// Restore CORE problem
	if( transform_SCEN_to_CORE(problem, change_and_restoration_point) != 0){
		printf("\nError while generating restoring CORE problem.\n");
		return(-1);
	}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nCORE problem restored.\n");
	
	// Free memory allocated on the for scenario differences
	free((*change_and_restoration_point).chgcoeff_col);
	free((*change_and_restoration_point).chgcoeff_row);
	free((*change_and_restoration_point).chgcoeff_coreval);
	free((*change_and_restoration_point).chgcoeff_scenval);
	free((*change_and_restoration_point).chgrhs_row);
	free((*change_and_restoration_point).chgrhs_coreval);
	free((*change_and_restoration_point).chgrhs_scenval);
	free((*change_and_restoration_point).chgbnd_col);
	free((*change_and_restoration_point).chgbnd_type);
	free((*change_and_restoration_point).chgbnd_coreval);
	free((*change_and_restoration_point).chgbnd_scenval);
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nRestoration point removed from memory.\n");
	
	// Restore output from Xpress to the value outside this function
	XPRSsetintcontrol(*problem, XPRS_OUTPUTLOG, xpress_output_log);
	
	// Return success indicator
	return(0);
}

/* Read data required by coordinator */
int read_SMPS_coordinator(const char *workdir, const char *SMPSrootname,
	struct string_buffer *firststagecolnames, double **update_scaling,
	struct string_buffer *scenarios, double **scenprobability,
	const struct options *opt)
{
	// Local variables that need to be freed at the end
	XPRSprob problem;
	struct string_buffer rownames, colnames, timestages, scenparent;
	int *rowststage, *colststage, *scenbranchtstage, *scenstartline, *scenendline, *firststagecols;
	double *fstg_lbnds, *fstg_ubnds;
	
	// Create problem
	XPRScreateprob(&problem);
	XPRSsetintcontrol(problem, XPRS_SCALING, 0);			// avoids automatic scaling when reading the problem
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		XPRSaddcbmessage(problem, Message, NULL, 0);	// defines callback function for printing log from Xpress (usefull only in Windows)
	#endif
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		XPRSsetintcontrol(problem, XPRS_OUTPUTLOG, 1);
	} else {
		XPRSsetintcontrol(problem, XPRS_OUTPUTLOG, 4);
	}
	
	// Load CORE problem to Xpress
	char mpsname[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nLoading CORE file to Xpress ...\n");
	strcpy(mpsname, workdir);
	strcat(strcat(strcat(mpsname, SYSTEM_SLASH), SMPSrootname), ".mps");
	XPRSreadprob(problem, mpsname, "");
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("Problem loaded in Xpress.\n");
	
	// Read problem row and column names
	int nl;
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nGathering row and column names from Xpress ...");
	XPRSgetintattrib(problem, XPRS_ROWS, &(rownames.num_elems));
	XPRSgetintattrib(problem, XPRS_COLS, &(colnames.num_elems));
	XPRSgetintattrib(problem, XPRS_NAMELENGTH, &nl);
	rownames.atom_len = uXPRS_NAMELENGTH_DEF*nl;
	colnames.atom_len = uXPRS_NAMELENGTH_DEF*nl;
	rownames.buffer = (char *) malloc(sizeof(char)*(rownames.atom_len+1)*rownames.num_elems);
	colnames.buffer = (char *) malloc(sizeof(char)*(colnames.atom_len+1)*colnames.num_elems);
	XPRSgetnames(problem, 1, rownames.buffer, 0, rownames.num_elems-1);
	XPRSgetnames(problem, 2, colnames.buffer, 0, colnames.num_elems-1);
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nRows and column names available to main program.\n");
	
	// Read TIME file
	char timename[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nParsing TIME file ...");
	strcpy(timename, workdir);
	strcat(strcat(strcat(timename, SYSTEM_SLASH), SMPSrootname), ".tim");
	timestages.atom_len = SMPS_MAXFIELDLEN;
	timestages.num_elems = 0;
	timestages.buffer = (char *) malloc(sizeof(char)*(timestages.atom_len+1)*SMPS_MAXSTAGES);
	rowststage = (int *) calloc(rownames.num_elems, sizeof(int));
	colststage = (int *) calloc(colnames.num_elems, sizeof(int));
	if( read_TIME_file(SMPSrootname, timename, &rownames, &colnames, &timestages, rowststage, colststage) ){
		printf("\nError while reading TIME file.\n");
		return(-1);
		}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\n%d time stages detected. All rows and cols correctly assigned to time stages.\n",
			timestages.num_elems);
	
	// Read general scenario info from STOCH file
	char stochname[200];
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\nPreprocessing STOCH file:\n(i) Getting number and name length of scenarios ...");
	strcpy(stochname, workdir);
	strcat(strcat(strcat(stochname, SYSTEM_SLASH), SMPSrootname), ".sto");
	
	// get scenario set size for creating buffers
	if( get_scenario_sizes(SMPSrootname, stochname, &((*scenarios).atom_len), &((*scenarios).num_elems)) ){
		printf("\nError while reading the character length and number of scenarios.\n");
		return(-1);
		}
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("\n%d scenarios with names up to %d characters long detected in STOCH file.\n",
			(*scenarios).num_elems, (*scenarios).atom_len);
	
	// global buffers
	(*scenarios).buffer = (char *) malloc(sizeof(char)*((*scenarios).atom_len+1)*(*scenarios).num_elems);
	*scenprobability = (double *) calloc((*scenarios).num_elems, sizeof(double));
	
	// local buffers
	scenparent.atom_len = (int) fmax((*scenarios).atom_len, strlen(SMPS_ROOTNODENAME));
	scenparent.num_elems = (*scenarios).num_elems;
	scenparent.buffer = (char *) malloc(sizeof(char)*(scenparent.atom_len+1)*scenparent.num_elems);
	scenbranchtstage = (int *) calloc((*scenarios).num_elems, sizeof(int));
	scenstartline = (int *) calloc((*scenarios).num_elems, sizeof(int));
	scenendline = (int *) calloc((*scenarios).num_elems, sizeof(int));
	
	// read scenario metadata into buffers
	if((*opt).verbosity >= VERB_WORKER_FUNCS)
		printf("(ii) Allocating space and reading scenario metadata ...");	
	if( get_scenario_metadata(SMPSrootname, stochname, &timestages,
		scenarios, &scenparent, scenbranchtstage, *scenprobability, scenstartline, scenendline) ){
		printf("\nError while reading scenarios' metadata.\n");
		return(-1);
		}
	
	// verify that the probabilities add up to 1
	{
		int i;
		double totalprob = 0;
		for(i = 0; i < (*scenarios).num_elems; i++){
			totalprob += *((*scenprobability) + i);
		}
		if( totalprob < 1.0){
			printf("\nWARNING: Scenario probalities are %e less than 1. They will be rescaled up to 1.\n", 1.0-totalprob);
			for(i = 0; i < (*scenarios).num_elems; i++){
				*((*scenprobability) + i) /= totalprob;
			}
		}
	}
	
	// print scenario data summary
	if((*opt).verbosity >= VERB_WORKER_FUNCS){
		int i, minchg, maxchg;
		float minprob, maxprob;
		minchg = *(scenendline + (*scenarios).num_elems - 1);
		maxchg = 0;
		minprob = 1;
		maxprob = 0;
		for(i = 0; i < (*scenarios).num_elems; i++){
			minchg = (int) fmin(minchg, (*(scenendline + i)) - (*(scenstartline + i)));
			maxchg = (int) fmax(maxchg, (*(scenendline + i)) - (*(scenstartline + i)));
			minprob = (float) fmin(minprob, *((*scenprobability) + i));
			maxprob = (float) fmax(maxprob, *((*scenprobability) + i));
		}
		printf("\nScenario metadata for %d scenarios read from STOCH file."
			"\nNumber of changes to CORE (min - max):\t%d - %d"
			"\nProbabilities (min - max):\t\t%.2f%% - %.2f%%\n",
			(*scenarios).num_elems, minchg, maxchg, 100*minprob, 100*maxprob);
	}
	
	// Determine first stage columns
	char firststagefilled[timestages.atom_len+1];
	add_spaces(SMPS_FIRSTSTAGENAME, timestages.atom_len+1, timestages.atom_len, firststagefilled);
	int firststageind = get_string_index(&timestages, &firststagefilled[0]);
	if(firststageind < 0){
		printf("\nFirst stage '" SMPS_FIRSTSTAGENAME "' not found among time stages.\n");
		return(-1);
	}
	firststagecols = (int*) malloc(sizeof(int)*colnames.num_elems);
	select_cols_and_rows(&problem, &rownames, &colnames, 0, 0, &timestages, rowststage, colststage,
		1, &firststageind, 0, NULL, &((*firststagecolnames).num_elems), firststagecols,
		NULL, NULL, NULL);
	
	// Write first stage cols global variable
	(*firststagecolnames).atom_len = colnames.atom_len;
	(*firststagecolnames).buffer = (char*) malloc(sizeof(char)*
		((*firststagecolnames).atom_len + 1)*(*firststagecolnames).num_elems);
	{
		int i;
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			insert_in_buffer(firststagecolnames, i,
				colnames.buffer+(colnames.atom_len+1) * *(firststagecols + i));
		}
	}
	
	// Get first stage cols bounds from CORE problem
	fstg_lbnds = (double*) malloc(sizeof(double)*(*firststagecolnames).num_elems);
	fstg_ubnds = (double*) malloc(sizeof(double)*(*firststagecolnames).num_elems);
	{
		int i;
		double *bndall = (double*) malloc(sizeof(double)*colnames.num_elems);
		XPRSgetlb(problem, bndall, 0, colnames.num_elems - 1);
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			*(fstg_lbnds + i) = *(bndall + *(firststagecols + i));
		}
		XPRSgetub(problem, bndall, 0, colnames.num_elems - 1);
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			*(fstg_ubnds + i) = *(bndall + *(firststagecols + i));
		}
		free(bndall);
	}
	
	// Gather outer bounds from STOCH file, scenario by scenario
	{
		int i;
		for(i = 0; i < (*scenarios).num_elems; i++){
			if( get_scenario_col_bounds(SMPSrootname, stochname,
				(*scenarios).buffer + ((*scenarios).atom_len + 1) * i,
				scenarios, &scenparent, scenstartline, scenendline,
				firststagecolnames, fstg_lbnds, fstg_ubnds) != 0 ){
				printf("\nProblem reading outer bounds for first stage variables at scenario %s.",
					(*scenarios).buffer + ((*scenarios).atom_len + 1) * i);
				return(-1);
			}
		}
	}
	
	// Check that all first stage variables are bounded
	{
		int i;
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			if( *(fstg_lbnds + i) <= XPRS_MINUSINFINITY || *(fstg_ubnds + i) >= XPRS_PLUSINFINITY ){
				printf("\nUnbounded column %s detected in first stage matrix.\n",
					(*firststagecolnames).buffer + ((*firststagecolnames).atom_len + 1) * i );
				return(-1);
			}
		}
	}
	
	/*
	// debug only
	// Print first stage variable bounds
	{
		FILE *file = fopen("bounds.txt", "w");
		fprintf(file, "\nVariable\tLB\tUB\n");
		int i;
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			fprintf(file, "%s\t%f\t%f\n",
				(*firststagecolnames).buffer + ((*firststagecolnames).atom_len + 1) * i,
				*(fstg_lbnds + i), *(fstg_ubnds + i));
		}
		fclose(file);
		return(-1);
	}
	// debug only
	*/
	
	// Compute scaling
	*update_scaling = (double*) malloc(sizeof(double)*(*firststagecolnames).num_elems);
	{
		int i;
		for(i = 0; i < (*firststagecolnames).num_elems; i++){
			if( (*(fstg_ubnds + i)) - (*(fstg_lbnds + i)) > DBL_EPSILON
				&& ((*(fstg_ubnds + i)) - (*(fstg_lbnds + i)) >= 1.0 + DBL_EPSILON
				|| (*(fstg_ubnds + i)) - (*(fstg_lbnds + i)) <= 1.0 - DBL_EPSILON) ){
				*((*update_scaling) + i) = 1.0/(*(fstg_ubnds + i)) - (*(fstg_lbnds + i));
			} else {
				*((*update_scaling) + i) = NO_UPDATE_SCALING;
			}
		}
	}
	
	// Delete problem and free heap
	XPRSdestroyprob(problem);
	free(rownames.buffer);
	free(colnames.buffer);
	free(timestages.buffer);
	free(scenparent.buffer);
	free(rowststage);
	free(colststage);
	free(scenbranchtstage);
	free(scenstartline);
	free(scenendline);
	free(firststagecols);
	free(fstg_lbnds);
	free(fstg_ubnds);
	
	// Return success indicator
	return(0);
}

/* Free memory from data needed by the coordinator process */
int free_heap_from_read_SMPS_coordinator(
	struct string_buffer *firststagecolnames, double **update_scaling,
	struct string_buffer *scenarios, double **scenprobability, const struct options *opt)
{
	// Free memory allocated on the HEAP
	free((*firststagecolnames).buffer);
	free((*scenarios).buffer);
	free(*scenprobability);
	free(*update_scaling);
	
	// Return success indicator
	return(0);
}

/* Function for generate MPI message types */
int create_message_types(const int numfirststagecols,
	MPI_Datatype *dualscenjobmsg, MPI_Datatype *dualnonscenjobmsg, MPI_Datatype *dualresultmsg,
	MPI_Datatype *primaljobmsg, MPI_Datatype *primalprojresultmsg, MPI_Datatype *primalresultmsg)
{
	// Elements used by the function call
	int count, *array_of_blocklengths;
	MPI_Aint *array_of_displacements;
	MPI_Datatype *array_of_types;
	MPI_Aint lb, extent_int, extent_double;
	MPI_Type_get_extent(MPI_INT, &lb, &extent_int);
	MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent_double);
	
	// Construct dual scenario job message type = (Int, Double[numfirststagecols])
	// Dual job: (ScenarioNum, x[])
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 1;
	*(array_of_blocklengths + 1) = numfirststagecols;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, dualscenjobmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI dual job message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Construct dual nonscenario job message type = (Int, Double[1+3*numfirststagecols])
	// Dual non scenario job: (NOSCENID, LBscen, msumx_delayed[], msumx_current[], ucenter[])
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 1;
	*(array_of_blocklengths + 1) = 1+3*numfirststagecols;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, dualnonscenjobmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI dual job message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Construct dual result message type: (Int[2], Double[3+numfirststagecols])
	// Dual result: (ScenarioNum, XPRS status, worktime, Obj (UB), LB, u[])
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 2;
	*(array_of_blocklengths + 1) = 3+numfirststagecols;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = 2*extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, dualresultmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI dual result message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Construct primal job message type = (Int[2], Double[numfirststagecols])
	// Primal job: (CandidateID, ScenarioNum, u[])
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 2;
	*(array_of_blocklengths + 1) = numfirststagecols;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = 2*extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, primaljobmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI primal job message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Construct primal projection result message type = (Int, Double[1+numfirststagecols])
	// Primal projection job result: (CandidateID, worktime, u[])
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 1;
	*(array_of_blocklengths + 1) = 1+numfirststagecols;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, primalprojresultmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI primal job message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Construct primal result message type: (Int[3], Double[3])
	// Primal result: (CandidateID, ScenarioNum, XPRS status, worktime, Obj, LB)
	count = 2;
	array_of_blocklengths = (int*) malloc(sizeof(int)*count);
	*array_of_blocklengths = 3;
	*(array_of_blocklengths + 1) = 3;
	array_of_displacements = (MPI_Aint*) malloc(sizeof(MPI_Aint)*count);
	*array_of_displacements = 0;
	*(array_of_displacements + 1) = 3*extent_int;
	array_of_types = (MPI_Datatype*) malloc(sizeof(MPI_Datatype)*count);
	*array_of_types = MPI_INT;
	*(array_of_types + 1) = MPI_DOUBLE;
	if( MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
		array_of_types, primalresultmsg) != MPI_SUCCESS ){
		printf("\nError creating MPI primal result message type\n.");
		return(-1);
	}
	free(array_of_blocklengths);
	free(array_of_displacements);
	free(array_of_types);
	
	// Commit MPI_Datatypes
	MPI_Type_commit(dualscenjobmsg);
	MPI_Type_commit(dualnonscenjobmsg);
	MPI_Type_commit(dualresultmsg);
	MPI_Type_commit(primaljobmsg);
	MPI_Type_commit(primalprojresultmsg);
	MPI_Type_commit(primalresultmsg);
	
	// Return success indicator
	return(0);
}

/* Function to free message types */
int free_message_types(
	MPI_Datatype *dualscenjobmsg, MPI_Datatype *dualnonscenjobmsg, MPI_Datatype *dualresultmsg,
	MPI_Datatype *primaljobmsg, MPI_Datatype *primalprojresultmsg, MPI_Datatype *primalresultmsg)
{
	// Free MPI_Datatypes created for communications
	MPI_Type_free(dualscenjobmsg);
	MPI_Type_free(dualnonscenjobmsg);
	MPI_Type_free(dualresultmsg);
	MPI_Type_free(primaljobmsg);
	MPI_Type_free(primalprojresultmsg);
	MPI_Type_free(primalresultmsg);
	
	// Return success indicator
	return(0);
}