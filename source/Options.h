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

#ifndef OPTIONS_INCLUDED
	
	// define OPTIONS_INCLUDED to avoid loading the header again
	#define OPTIONS_INCLUDED
	
	/* INCLUDE STATEMENTS */
	
	// Standard libraries
	#include <stdio.h>			// provides functions for reading and printing
	#include <stdlib.h>			// provides several functions for managing memory, files and interact with system
	#include <string.h>			// managing strings
	
	/* CONSTANTS */
	
	// Configuration file constants
	#define CONFIG_NAMECOL			4
	#define CONFIG_VALUECOL			24
	#define CONFIG_MAXLINELEN		1000
	#define CONFIG_MAXNAMELEN		18
	#define CONFIG_MAXNAMELENSTR		"18"
	#define CONFIG_VERBOSITY		"VERBOSITY"
	#define CONFIG_MAXITER			"MAXITER"
	#define CONFIG_MAXTIME			"MAXTIME"
	#define CONFIG_DUAL_SHARE_MIDITER	"DSHARE_MIDITER"		// mid point for dual share change
	#define CONFIG_DUAL_SHARE_START		"DSHARE_START"			// start dual share
	#define CONFIG_DUAL_SHARE_END		"DSHARE_END"			// end dual share
	#define CONFIG_NONSCEN_PERIOD		"NONSCEN_PERIOD"		// evaluation of f0 every how many dual updates (in pu of number of scenarios)
	#define CONFIG_NONSCEN_CENTER		"NONSCEN_CENTER"		// center type for smoothing: fixed center during initializations, moving center with the incumbent or moving center with subgradients (ADMM like)
	#define CONFIG_DUAL_STEPSIZE		"DUAL_STEPSIZE"
	#define CONFIG_STEPSIZE_P		"DUAL_STEPSIZE_P"
	#define CONFIG_STEPSIZE_R		"DUAL_STEPSIZE_R"
	#define CONFIG_STEPSIZE_Q		"DUAL_STEPSIZE_Q"
	#define CONFIG_STEPSIZE_XI		"DUAL_STEPSIZE_XI"
	#define CONFIG_STEPSIZE_SGM		"DUAL_STEPSIZE_SGM"
	#define CONFIG_GNORM_FILTER		"DUAL_GNORM_FILTER"		// Filtering the norm of the subgradient for Polyak stepsie (zero corresponds to Polyak). See Camerini, Fratta, and Maffioli 75.
	#define CONFIG_PRIMAL_RECOVERY		"PRIMAL_RECOVERY"
	#define CONFIG_PRIMAL_IS_SSIZE		"PRIMAL_IS_SSIZE"		// IS sample size for primal recovery
	#define CONFIG_INIT_TYPE		"INIT_TYPE"				// Initialization type: linear relaxation or period relaxation
	#define CONFIG_RELDUALTOL		"DUAL_GAP_TOL"
	#define CONFIG_SMOOTHING_PARAM		"SMOOTHING_PARAM"
	#define CONFIG_PROJ_ROUND_MIDPOINT	"PROJ_ROUND_MIDP"		// mid point for rounding when projecting infeasible first stage solutions (lower is more conservative)
	#define CONFIG_PERIOD_AGG_POLICY	"PERIOD_AGG_POL"		// policy for aggregating "multi-period" first stage variables in period relaxation
	#define CONFIG_SCEN_LP_ALG		"SCEN_LP_ALG"
	#define CONFIG_PERIOD_LP_ALG		"PERIOD_LP_ALG"
	
	// Options constants
	#define VERB_NONE			0
	#define VERB_COORDINATOR_ONLY		1
	#define VERB_COORDINATOR_FUNCS		2
	#define VERB_WORKER			3
	#define VERB_WORKER_FUNCS		4
	#define CENTER_FIXED			10
	#define CENTER_INCUMBENT		11
	#define CENTER_AVERAGE_SUBGRADIENT	12
	#define DUAL_STEPSIZE_CONSTANT		20		// p
	#define DUAL_STEPSIZE_DIMINISHING	21		// p/(1+r/numscenarios*k)^q
	#define DUAL_STEPSIZE_POLYAK		22		// p/(1+r/numscenarios*k)^q*min(xi*max(abs(f(x)),abs(fup)), fup - f(x))/max(sigma, |g|^2_2)
	#define PRIMAL_RECOVERY_FIFO		30		// getting solutions from scenario subproblems as they appear -> requires to keep all solutions
	#define PRIMAL_RECOVERY_RND		31		// getting the solution of one scenario at random from the pile -> requires to keep all solutions
	#define PRIMAL_RECOVERY_LIFO		32		// getting the solution of the lastly solved scenario
	#define PRIMAL_RECOVERY_IS		33		// Importance Sampling recombination heuristic
	#define INIT_TYPE_LINEAR		40		// Use linear relaxation without delayed rows for initialization
	#define INIT_TYPE_PERIOD		41		// Use period relaxation for initialization
	#define PERIOD_POLICY_MIN		50		// Aggregate "multi-period" first stage variables using min 
	#define PERIOD_POLICY_MEAN		51		// Aggregate "multi-period" first stage variables using mean
	#define PERIOD_POLICY_MAX		52		// Aggregate "multi-period" first stage variables using max
	
	/* STRUCTS */
	typedef struct xpress_options {
		char lp_method[2];
		int num_intoptions, num_dbloptions;
		int *intoptions, *dbloptions;
		int *intoptions_value;
		double *dbloptions_value;
	} xpress_options;
	typedef struct options {
		int verbosity, dual_stepsize, primal_recovery;			// caterigorical
		int max_iterations, init_type, period_aggregation_policy;
		double max_time, dual_share_midpoint, dual_share_start, dual_share_end,
			nonscen_evaluation_period, nonscen_center_type;
		double dual_stepsize_p, dual_stepsize_r, dual_stepsize_q, dual_stepsize_xi, dual_stepsize_sigma,
			dual_gnorm2_filter;
		double primal_IS_ssize, rel_dual_gap_tol, dual_f0_mu, proj_round_midpoint;
		xpress_options scenprob_options, periodprob_options;
	} options;
	
	/* FUNCTIONS AND PROCEDURES */
	extern int read_options(const char*, struct options*);
	extern int free_options(struct options*);
	
#endif
