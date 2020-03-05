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
 * Read options from configuration file
*/

/* Header files */
#include "Options.h"
#include "LowLevelFunc.h"

/* Prototypes */
int read_xpress_intoptions(FILE *file, const int numoptions, int *intoptions, int *intoptions_value);
int read_xpress_dbloptions(FILE *file, const int numoptions, int *dbloptions, double *dbloptions_value);

/* Read options from configuration file */
int read_options(const char *optionsfile, struct options *globalopt)
{
	// variables
	char line[CONFIG_MAXLINELEN], varname[CONFIG_MAXNAMELEN];
	int lncnt = 0;
	
	// open configuration file for reading
	FILE *file = fopen(optionsfile, "r");
	if (file == NULL) {
		printf("\nConfiguration file does not exist.\n");
		return(-1);									// error: file does not exist
	}
	
	// skip lines until finding GENERAL
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, "GENERAL") != 0){
		lncnt++;
	}
	lncnt++;
	
	// reading verbosity level
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_VERBOSITY) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).verbosity)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading verbosity level from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading max iterations
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_MAXITER) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).max_iterations)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading max iterations from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading max time
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_MAXTIME) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).max_time)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading max time from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading mid point for dual share change
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_DUAL_SHARE_MIDITER) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_share_midpoint)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading mid point for dual share change from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual share start value
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_DUAL_SHARE_START) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_share_start)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual share start value from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual share end value
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_DUAL_SHARE_END) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_share_end)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual share end value from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading non-scenario subproblem execution period
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_NONSCEN_PERIOD) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).nonscen_evaluation_period)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading non scenario subproblem evaluation period from %s."
		"\nLine: %s\n", optionsfile, line);
		return(-1);
	}
	
	// reading non-scenario subproblem execution period
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_NONSCEN_CENTER) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).nonscen_center_type)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading non scenario subproblem center type from %s."
		"\nLine: %s\n", optionsfile, line);
		return(-1);
	}
	
	// reading dual stepsize type
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_DUAL_STEPSIZE) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).dual_stepsize)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize type from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual stepsize p parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_STEPSIZE_P) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_stepsize_p)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize parameter p from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual stepsize r parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_STEPSIZE_R) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_stepsize_r)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize parameter r from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual stepsize q parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_STEPSIZE_Q) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_stepsize_q)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize parameter q from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual stepsize xi parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_STEPSIZE_XI) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_stepsize_xi)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize parameter xi from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading dual stepsize sigma parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_STEPSIZE_SGM) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_stepsize_sigma)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading dual stepsize parameter sigma from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading gnorm^2 filter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_GNORM_FILTER) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_gnorm2_filter)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading subgradient filer parameter from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading primal recovery type
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_PRIMAL_RECOVERY) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).primal_recovery)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading primal recovery type from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading sample size for IS primal recovery
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_PRIMAL_IS_SSIZE) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).primal_IS_ssize)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while sample size for IS primal recovery from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading initialization type
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_INIT_TYPE) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).init_type)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading initialization type from %s.\n", optionsfile);
		return(-1);
	}
	
	// reading duality gap tolerance
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_RELDUALTOL) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).rel_dual_gap_tol)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading duality gap tolerance from %s.\n", optionsfile);
		return(-1);
	}
	
	// smoothing parameter
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_SMOOTHING_PARAM) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).dual_f0_mu)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading duality gap tolerance from %s.\n", optionsfile);
		return(-1);
	}
	
	// rounding midpoint
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_PROJ_ROUND_MIDPOINT) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%lf", &((*globalopt).proj_round_midpoint)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading rounding mid point from %s.\n", optionsfile);
		return(-1);
	}
	
	// policy for aggregating "multi-period" first stage variables in period relaxation
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_PERIOD_AGG_POLICY) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%d", &((*globalopt).period_aggregation_policy)) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading period aggregation policy from %s.\n", optionsfile);
		return(-1);
	}
	
	// scenario subproblem algorithm
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_SCEN_LP_ALG) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%1s", &(((*globalopt).scenprob_options).lp_method)[0]) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading scenario LP algorithm from %s.\n", optionsfile);
		return(-1);
	}
	
	// period subproblem algorithm
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line+CONFIG_NAMECOL, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, CONFIG_PERIOD_LP_ALG) == 0
		&& sscanf(line+CONFIG_VALUECOL, "%1s", &(((*globalopt).periodprob_options).lp_method)[0]) == 1){
		lncnt++;
	} else {
		printf("\nProblem found while reading period LP algorithm from %s.\n", optionsfile);
		return(-1);
	}
	
	// check consistency of double valued options
	if( ((*globalopt).dual_share_midpoint <= 0.0 || (*globalopt).dual_share_midpoint >= 1.0 )
		|| ((*globalopt).dual_share_start <= 0.0 || (*globalopt).dual_share_start >= 1.0 )
		|| ((*globalopt).dual_share_end <= 0.0 || (*globalopt).dual_share_end >= 1.0 )
		|| ((*globalopt).nonscen_evaluation_period <= 0.0 || (*globalopt).nonscen_evaluation_period >= 1.0 )
		|| (*globalopt).dual_stepsize_p < 0
		|| (*globalopt).dual_stepsize_r < 0
		|| ((*globalopt).dual_stepsize_q <= 0.5 || (*globalopt).dual_stepsize_q > 1.0 )
		|| ((*globalopt).dual_stepsize_xi <= 0.0 || (*globalopt).dual_stepsize_xi >= 1.0 )
		|| ((*globalopt).primal_IS_ssize <= 0.0 || (*globalopt).primal_IS_ssize >= 1.0 )
		|| (*globalopt).dual_f0_mu <= 0
		|| ((*globalopt).proj_round_midpoint <= 0.0 || (*globalopt).proj_round_midpoint >= 1.0 ) ){
		printf("\nInconsisteng parameters were read from %s."
			"\nDual share midpoint, start and end values must be within (0,1)."
			"\nNon-scenario evaluation period," CONFIG_NONSCEN_PERIOD ", must be within (0,1)."
			"\nDual stepsize parameters must all be non-negative, 0.5 < q <= 1 and 0 < xi <= 1."
			"\nSample size for IS primal recovery," CONFIG_PRIMAL_IS_SSIZE ", must be within (0,1)."
			"\nSmoothing parameter," CONFIG_SMOOTHING_PARAM ", must be positive."
			"\nRounding mid pont," CONFIG_PROJ_ROUND_MIDPOINT ", must be within (0,1).\n", optionsfile);
		return(-1);
	}
	
	// get lenght of each xpress parameter section
	int xpress_params_ix[5];
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, "XPRESS_SCEN_INT") != 0){
		lncnt++;
	}
	xpress_params_ix[0] = lncnt;
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, "XPRESS_SCEN_DBL") != 0){
		lncnt++;
	}
	xpress_params_ix[1] = lncnt;
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, "XPRESS_PERIOD_INT") != 0){
		lncnt++;
	}
	xpress_params_ix[2] = lncnt;
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" CONFIG_MAXNAMELENSTR "s", varname) == 1
		&& strcmp(varname, "XPRESS_PERIOD_DBL") != 0){
		lncnt++;
	}
	xpress_params_ix[3] = lncnt;
	while(fgets(line,sizeof line,file) !=  NULL){
		lncnt++;
	}
	xpress_params_ix[4] = lncnt;
	
	// close configuration file
	fclose(file);
	
	// allocate buffers to contain xpress configuration for scenario subproblems
	(*globalopt).scenprob_options.num_intoptions = xpress_params_ix[1] - xpress_params_ix[0];
	(*globalopt).scenprob_options.intoptions = 
		(int*) malloc( sizeof(int) * (*globalopt).scenprob_options.num_intoptions );
	(*globalopt).scenprob_options.intoptions_value = 
		(int*) malloc( sizeof(int) * (*globalopt).scenprob_options.num_intoptions );
	(*globalopt).scenprob_options.num_dbloptions = xpress_params_ix[2] - xpress_params_ix[1];
	(*globalopt).scenprob_options.dbloptions = 
		(int*) malloc( sizeof(int) * (*globalopt).scenprob_options.num_dbloptions );
	(*globalopt).scenprob_options.dbloptions_value = 
		(double*) malloc( sizeof(double) * (*globalopt).scenprob_options.num_dbloptions );
	
	// allocate buffers to contain xpress configuration for period subproblems
	(*globalopt).periodprob_options.num_intoptions = xpress_params_ix[3] - xpress_params_ix[2];
	(*globalopt).periodprob_options.intoptions = 
		(int*) malloc( sizeof(int) * (*globalopt).periodprob_options.num_intoptions );
	(*globalopt).periodprob_options.intoptions_value = 
		(int*) malloc( sizeof(int) * (*globalopt).periodprob_options.num_intoptions );
	(*globalopt).periodprob_options.num_dbloptions = xpress_params_ix[4] - xpress_params_ix[3];
	(*globalopt).periodprob_options.dbloptions = 
		(int*) malloc( sizeof(int) * (*globalopt).periodprob_options.num_dbloptions );
	(*globalopt).periodprob_options.dbloptions_value = 
		(double*) malloc( sizeof(double) * (*globalopt).periodprob_options.num_dbloptions );
	
	// re-open file and skip lines up to XPRESS_SCEN_INT
	file = fopen(optionsfile, "r");
	flineskip(file, xpress_params_ix[0], CONFIG_MAXLINELEN);
	
	// collect scenario subproblem Xpress parameters
	if( read_xpress_intoptions(file, (*globalopt).scenprob_options.num_intoptions,
		(*globalopt).scenprob_options.intoptions, (*globalopt).scenprob_options.intoptions_value) != 0 ||
		read_xpress_dbloptions(file, (*globalopt).scenprob_options.num_dbloptions,
		(*globalopt).scenprob_options.dbloptions, (*globalopt).scenprob_options.dbloptions_value) != 0 ){
		printf("\nError while reading Xpress parameters for scenario subproblems.\n");
		return(-1);
	}
	
	// collect period subproblem Xpress parameters
	if( read_xpress_intoptions(file, (*globalopt).periodprob_options.num_intoptions,
		(*globalopt).periodprob_options.intoptions, (*globalopt).periodprob_options.intoptions_value) != 0 ||
		read_xpress_dbloptions(file, (*globalopt).periodprob_options.num_dbloptions,
		(*globalopt).periodprob_options.dbloptions, (*globalopt).periodprob_options.dbloptions_value) != 0 ){
		printf("\nError while reading Xpress parameters for period subproblems.\n");
		return(-1);
	}
	
	// return success indicator
	return(0);
}

/* Release space used to store Xpress options */
int free_options(struct options *globalopt)
{
	// free malloc calls made by read_options
	free((*globalopt).scenprob_options.intoptions);
	free((*globalopt).scenprob_options.intoptions_value);
	free((*globalopt).scenprob_options.dbloptions);
	free((*globalopt).scenprob_options.dbloptions_value);
	free((*globalopt).periodprob_options.intoptions);
	free((*globalopt).periodprob_options.intoptions_value);
	free((*globalopt).periodprob_options.dbloptions);
	free((*globalopt).periodprob_options.dbloptions_value);
	
	// return success indicator
	return(0);
}

/* Read integer xpress control params */
int read_xpress_intoptions(FILE *file, const int numoptions, int *intoptions, int *intoptions_value)
{
	// skip section header
	flineskip(file, 1, CONFIG_MAXLINELEN);
	
	// parse lines one by one
	char line[CONFIG_MAXLINELEN];
	int i;
	for(i = 0; i < numoptions; i++){
		if( !( fgets(line,sizeof line,file) !=  NULL
			&& sscanf(line+CONFIG_NAMECOL, "%d", intoptions + i) == 1
			&& sscanf(line+CONFIG_VALUECOL, "%d", intoptions_value + i) == 1 ) ){
			printf("\nProgram failed at collecting integer control parameters for Xpress.\n");
			printf("%s", line);
			return(-1);
		}
	}
	
	// return success indicator
	return(0);
}

/* Read double xpress control params */
int read_xpress_dbloptions(FILE *file, const int numoptions, int *dbloptions, double *dbloptions_value)
{
	// skip section header
	flineskip(file, 1, CONFIG_MAXLINELEN);
	
	// parse lines one by one
	char line[CONFIG_MAXLINELEN];
	int i;
	for(i = 0; i < numoptions; i++){
		if( !( fgets(line,sizeof line,file) !=  NULL
			&& sscanf(line+CONFIG_NAMECOL, "%d", dbloptions + i) == 1
			&& sscanf(line+CONFIG_VALUECOL, "%lf", dbloptions_value + i) == 1 ) ){
			printf("\nProgram failed at collecting double control parameters for Xpress.\n");
			return(-1);
		}
	}
	
	// return success indicator
	return(0);
}

/* Configuration file structure
 * Fixed column positions: names 5 and values 25
 * Order of sections must be respected. Order of params within general section must also be respected.
GENERAL # this will be ignored
    OPTNAME            VALUE # this will be ignored
XPRESS_SCEN_INT # this will be ignored
    XPRSID             VALUE  # this will be ignored
XPRESS_SCEN_DBL # this will be ignored
    XPRSID             VALUE  # this will be ignored
XPRESS_PERIOD_INT # this will be ignored
    XPRSID             VALUE   # this will be ignored
XPRESS_PERIOD_DBL # this will be ignored
    XPRSID             VALUE   # this will be ignored
*/