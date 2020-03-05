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

#ifndef SMPS_XPRESS_INCLUDED
	
	// define SMPS_XPRESS_INCLUDED to avoid loading the header again
	#define SMPS_XPRESS_INCLUDED
	
	/* INCLUDE STATEMENTS */
	
	// Standard libraries
	#include <stdio.h>			// provides functions for reading and printing
	#include <stdlib.h>			// provides several functions for managing memory, files and interact with system
	#include <stddef.h>			// definitions (apparently needed to define size_t)
	#include <sys/stat.h>		// allows to manage files and directories
	#include <string.h>			// managing strings
	
	// Xpress API
	#include <xprs.h>			// provides functions for managing Xpress from C, C++
	
	// Custom header files
	#include "LowLevelFunc.h"
	#include "Options.h"
	
	// Xpress API
	//#include <xprs.h>			// provides functions for managing Xpress from C, C++
	
	/* CONSTANTS */
	
	// SMPS constants
	#define SMPS_MAXLINELEN				1000	// overly conservative to avoid surprises
	#define SMPS_MAXLINELENSTR			"1000"	// C sucks!
	#define SMPS_CODELEN				2		// 2 characters
	#define SMPS_CODELENSTR				"2"		// C sucks!
	#define SMPS_NOCODECHAR				"  "
	#define SMPS_MAXFIELDLEN			8		// 8 characters
	#define SMPS_MAXFIELDLENSTR			"8"		// C sucks!
	#define SMPS_MAXSTAGES				100		// overly conservative to avoid surprises
	#define SMPS_CODECOL				1
	#define SMPS_FIRSTNAMECOL			4
	#define SMPS_SECONDNAMECOL			14
	#define SMPS_FIRSTNUMCOL			24
	#define SMPS_THIRDNAMECOL			39
	#define SMPS_SECONDNUMCOL			49
	#define SMPS_ROOTNODENAME			"'ROOT'"
	#define SMPS_RHSNAME				"rhs"
	#define SMPS_RHSNAME_ALT			"RHS00001"
	#define SMPS_BNDNAME				"bnd"
	#define SMPS_BNDNAME_ALT			"BND00001"
	#define SMPS_FIRSTSTAGENAME			"T0"
	#define SMPS_SCENARIO_CODE			" SC "
	#define SMPS_ENDDATA				"ENDATA"
	
	// Constants used in calls to Xpress
	#define uXPRS_NAMELENGTH_DEF			8		// see definition of NAMELENGTH in Xpress manual
	
	// Other constants
	#define DBL_EPSILON				1E-8	 	// epsilon for double comparison
	
	// 7z executable directory
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		#define WIN_7ZIP_EXECUTABLE		"C:\\PROGRA~1\\7-Zip\\7z.exe -y"
	#endif
	
	// System specific syntax
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		#define SYSTEM_SLASH			"\\"
	#else
		#define SYSTEM_SLASH			"/"
	#endif
	
	/* STRUCTS */
	typedef struct scen_differences
	{
		int numchgcoeff, numchgrhs, numchgbnd;
		int *chgcoeff_col, *chgcoeff_row, *chgrhs_row, *chgbnd_col;
		char *chgbnd_type;
		double *chgcoeff_coreval, *chgrhs_coreval, *chgbnd_coreval;
		double *chgcoeff_scenval, *chgrhs_scenval, *chgbnd_scenval;
	} scen_differences;
	
	/* FUNCTIONS AND PROCEDURES */
	
	// Functions in SMPS_Prepare
	extern int prepare_files(const char*, const char*, const char*, const char*,
		const struct options*);
	
	// Functions in SMPS_Xpress_Read
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		extern void XPRS_CC Message(XPRSprob, void*, const char*, int, int);
	#endif
	extern int read_TIME_file(const char*, const char*,
		const struct string_buffer*, const struct string_buffer*,
		struct string_buffer*, int*, int*);
	extern int get_scenario_sizes(const char*, const char*, int*, int*);
	extern int get_scenario_metadata(const char*, const char*, const struct string_buffer*,
		struct string_buffer*, struct string_buffer*, int*, double*, int*, int*);
	extern int read_scen_differences(const char*, const char*, const char*,
		const struct string_buffer*, const struct string_buffer*, const int*,
		const struct string_buffer*, const struct string_buffer*,
		struct scen_differences*);
	extern int get_scenario_col_bounds(const char*, const char*, const char*,
		const struct string_buffer*, const struct string_buffer*,
		const int*, const int*, const struct string_buffer*, double*, double*);
	extern int read_delayed_rows(const char*, const struct string_buffer*, int*, int*);
	
	// Functions in SMPS_Xpress_ScenSubprob
	extern int transform_CORE_to_SCEN(XPRSprob*, struct scen_differences*);
	extern int transform_SCEN_to_CORE(XPRSprob*, const struct scen_differences*);
	extern int load_scenario_multipliers(XPRSprob*, const int, const int*,
		const double*, const double*, const double*, double*);
	extern int unload_scenario_multipliers(XPRSprob*, const int, const int*, const double*);
	extern int fix_first_stage_vars(XPRSprob*, const int, const int*,
		double*, double*, double*);
	extern int unfix_first_stage_vars(XPRSprob*, const int, const int*,
		const double*, const double*);
	extern int solve_sequential_period_relaxation(const XPRSprob*,
		const struct string_buffer*, const struct string_buffer*,
		const int, const int, const struct string_buffer*,
		const int*, const int*, const int, const int, const int*,
		double*, double*, const struct options*);
	
	// Functions in SMPS_Xpress_StageSubprob
	extern int select_cols_and_rows(const XPRSprob*, const struct string_buffer*,
		const struct string_buffer*, const int, const int, const struct string_buffer*,
		const int*, const int*, const int, int*, const int, int*, int*, int*,
		int*, int*, int*);
	extern int get_firststagecols_multiplicity(const XPRSprob*,
		const struct string_buffer*, const struct string_buffer*, const struct string_buffer*,
		const int*, const int*, const int, const int, const int*, int*);
	extern int load_params_period_subprob(const XPRSprob*, const char*,
		const struct string_buffer*, const struct string_buffer*,
		const int, const int, const struct string_buffer*,
		const int*, const int*, const int, const int, const int, const int*, const int*,
		XPRSprob*, int*, int*, int*);
	extern int load_first_stage_matrix(const XPRSprob*,
		const struct string_buffer*, const int*, const int, const int, const int*,
		XPRSprob*, const char*);
	extern int load_f0_quad_objective(const int, const struct string_buffer*,
		const int, const int*, const double*, const double*, const double*,
		const double, const double*, XPRSprob*, double*);
	extern int load_f0_lin_objective(const int, const double*, const double*,
		const double*, XPRSprob*);
	extern int clean_f0_objective(XPRSprob*);
	
#endif
