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
 * Functions for reading SMPS files into Xpress, using the C API
 */

/* SMPS Xpress algorithm header */
#include "SMPS_Xpress.h"

/* Function to produce output from Xpress in Windows */
#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
	void XPRS_CC Message(XPRSprob my_prob, void* object, const char *msg, int len, int msgtype)
	{
		switch(msgtype)
		{
			case 4: /* error */
			case 3: /* warning */
			case 2: /* not used */
			case 1: /* information */
			printf("%s\n", msg);
			break;
			default: /* exiting - buffers need flushing */
			fflush(stdout);
			break;
		}
	}
#endif

/* Function to read TIME file information */
int read_TIME_file(const char *problemname, const char *filename,
	const struct string_buffer *rows, const struct string_buffer *cols,			// input variables
	struct string_buffer *tstages, int *time_stage_row, int *time_stage_col)	// output variables
{
	char line[SMPS_MAXLINELEN+1];
	char field1[SMPS_MAXFIELDLEN+1], field2[SMPS_MAXFIELDLEN+1], longfield[SMPS_MAXLINELEN+1];
	char TIMEftype[SMPS_MAXFIELDLEN+1];
	struct string_buffer TIMEnames;
	int *TIMEindices, *TIMEsort, *COREsort;
	int i;

	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("\nTIME file does not exist.\n");
		return(-1);									// error: file does not exist
	}

	// check that first line correctly declares the file as a TIME file
	if(fgets(line, sizeof line, file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field1) == 1
		&& strcmp(field1, "TIME") == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(longfield, problemname) == 0){
		;
	} else {
		printf("\nWrong TIME file header.\n");
		return(-1);									// error: wrong file header
	}

	// read TIME file type: IMPLICIT or EXPLICIT
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field1) == 1
		&& strcmp(field1, "PERIODS") == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", TIMEftype) == 1
		&& (strcmp(TIMEftype, "IMPLICIT") == 0 || strcmp(TIMEftype, "EXPLICIT") == 0)){
		;
	} else {
		printf("\nUnknown PERIOD type declaration: %s\n", TIMEftype);
		return(-1);									// error: unknown period declaration
	}

	// read body of TIME file according to fixed width SMPS format [Gassmann2007]
	if(strcmp(TIMEftype, "EXPLICIT") == 0){
		// read rows header
		if(fgets(line, sizeof line, file) !=  NULL
			&& strcmp(line, "ROWS\n") == 0){
			;
		} else {
			printf("\nWrong header on line 3 of TIME file. Expected 'ROWS', found '%s'.\n", field1);
			return(-1);								// error: wrong time file
		}
		// read rows section into buffers
		TIMEnames.num_elems = (*rows).num_elems-1;						// OBJ does not appears in TIME file
		TIMEnames.atom_len = SMPS_MAXFIELDLEN;
		TIMEnames.buffer = (char *) malloc(sizeof(char)*(TIMEnames.atom_len+1)*TIMEnames.num_elems);
		TIMEindices = (int *) malloc(TIMEnames.num_elems*sizeof(int));
		i = 0;
		while(fgets(line,sizeof line,file) !=  NULL
			&& strcmp(line, "COLUMNS\n") != 0
			&& i < TIMEnames.num_elems){
			// extracting row name
			sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field1);
			if( insert_in_buffer(&TIMEnames, i, &field1[0]) != 0 ){
				printf("\nProblem detected while reading row names from TIME file.\n");
				return(-1);
			}
			// extracting period name
			sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2);
			*(TIMEindices + i) = get_index_within_buffer(tstages, &field2[0]);
			if(*(TIMEindices + i) < 0){
				printf("\nProblem detected while reading row time stages from TIME file.\n");
				return(-1);
			}
			i++;
		}
		// check whether we have time stage information for each row
		if(i < TIMEnames.num_elems && strcmp(line, "COLUMNS\n") == 0){
			printf("\nMissing time information for %d rows.\n", TIMEnames.num_elems - i);
			return(-1);
		}
		if(i==TIMEnames.num_elems && strcmp(line, "COLUMNS\n") != 0){
			printf("\nTIME file contains information for rows absent from CORE file.\n");
			return(-1);
		}
		// sort rows using quick sort from stdlib
		TIMEsort = (int *) malloc(TIMEnames.num_elems*sizeof(int));
		for(i = 0; i < TIMEnames.num_elems; i++) *(TIMEsort + i) = i;
		sort_r(TIMEsort, TIMEnames.num_elems, sizeof(int), &compare_atoms_buffer, &TIMEnames);
		COREsort = (int *) malloc(((*rows).num_elems - 1)*sizeof(int));
		for(i = 0; i < (*rows).num_elems - 1; i++) *(COREsort + i) = i+1;		// objective is the first row in the CORE file!
		sort_r(COREsort, (*rows).num_elems - 1, sizeof(int), &compare_atoms_buffer, (void*) rows);
		// assign rows to time stages
		*time_stage_row = -1;
		for(i = 0; i < (*rows).num_elems-1; i++){
			if( strncmp((*rows).buffer + ((*rows).atom_len+1) * *(COREsort + i),
				TIMEnames.buffer + (TIMEnames.atom_len+1) * *(TIMEsort + i),
				(int) fmin((*rows).atom_len, TIMEnames.atom_len)) != 0 ){
				printf("\nMissmatch between row names in CORE and TIME files (%d correct matchs):", i);
				printf("\n\t(%d)'%s' <> (%d)'%s'\t(%d chars)\n",
					*(COREsort + i),
					(*rows).buffer + ((*rows).atom_len+1) * *(COREsort + i),
					*(TIMEsort + i),
					TIMEnames.buffer + (TIMEnames.atom_len+1) * *(TIMEsort + i),
					(int) fmin((*rows).atom_len, TIMEnames.atom_len));
				return(-1);
			}
			*(time_stage_row + *(COREsort + i)) = *(TIMEindices + *(TIMEsort + i));
		}
		// release memory used to host row names and times
		free(TIMEindices);
		free(TIMEnames.buffer);
		free(COREsort);
		free(TIMEsort);
		// read time stage of each column
		TIMEnames.num_elems = (*cols).num_elems;
		TIMEnames.atom_len = SMPS_MAXFIELDLEN;
		TIMEnames.buffer = (char *) malloc(sizeof(char)*(TIMEnames.atom_len+1)*TIMEnames.num_elems);
		TIMEindices = (int *) malloc(TIMEnames.num_elems*sizeof(int));
		i = 0;
		while(fgets(line,sizeof line,file) !=  NULL
			&& strcmp(line, "ENDATA\n") != 0
			&& i < TIMEnames.num_elems){
			// extracting col name
			sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field1);
			if( insert_in_buffer(&TIMEnames, i, &field1[0]) != 0 ){
				printf("\nProblem detected while reading column names from TIME file.\n");
				return(-1);
			}
			// extracting period name
			sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2);
			*(TIMEindices + i) = get_index_within_buffer(tstages, &field2[0]);
			if(*(TIMEindices + i) < 0){
				printf("\nProblem detected while reading column time stages from TIME file.\n");
				return(-1);
			}
			i++;
		}
		// check whether we have time stage information for each column
		if(i < TIMEnames.num_elems && strcmp(line, "ENDATA\n") == 0){
			printf("\nMissing time information for %d columns.\n", TIMEnames.num_elems - i);
			return(-1);
		}
		if(i==TIMEnames.num_elems && strcmp(line, "ENDATA\n") != 0){
			printf("\nTIME file contains information for columns absent from CORE file.\n");
			return(-1);
		}
		// sort columns using quick sort from stdlib
		TIMEsort = (int *) malloc(TIMEnames.num_elems*sizeof(int));
		for(i = 0; i < TIMEnames.num_elems; i++) *(TIMEsort + i) = i;
		sort_r(TIMEsort, TIMEnames.num_elems, sizeof(int), &compare_atoms_buffer, &TIMEnames);
		COREsort = (int *) malloc((*cols).num_elems*sizeof(int));
		for(i = 0; i < (*cols).num_elems; i++) *(COREsort + i) = i;
		sort_r(COREsort, (*cols).num_elems, sizeof(int), &compare_atoms_buffer, (void*) cols);
		// assign columns to time stages
		for(i = 0; i < (*cols).num_elems; i++){
			if( strncmp((*cols).buffer + ((*cols).atom_len+1) * *(COREsort + i),
				TIMEnames.buffer + (TIMEnames.atom_len+1) * *(TIMEsort + i),
				(int) fmin((*cols).atom_len, TIMEnames.atom_len)) != 0 ){
				printf("\nMissmatch between column names at CORE and TIME files:");
				printf("\n\t(%d)'%s' <> (%d)'%s'\t(%d chars)\n",
					*(COREsort + i),
					(*cols).buffer + ((*cols).atom_len+1) * *(COREsort + i),
					*(TIMEsort + i),
					TIMEnames.buffer + (TIMEnames.atom_len+1) * *(TIMEsort + i),
					(int) fmin((*cols).atom_len, TIMEnames.atom_len));
				return(-1);
			}
			*(time_stage_col + *(COREsort + i)) = *(TIMEindices + *(TIMEsort + i));
		}
		// release memory used to host column names and times
		free(TIMEindices);
		free(TIMEnames.buffer);
		free(COREsort);
		free(TIMEsort);
	} else {
		printf("\nTIME file type %s description not currently implemented.\n", TIMEftype);
		return(-1);								// error: pending implementation of implicit
	}
	// close connection with file
	fclose(file);

	// Return success indicator
	return(0);
}

/* Function to get the name length and the number of scenarios */
int get_scenario_sizes(const char *problemname, const char *filename, int *scenarionamelen, int *numscenarios)
{
	char line[SMPS_MAXLINELEN+1];
	char field[SMPS_MAXFIELDLEN+1], longfield[SMPS_MAXLINELEN+1];

	// open connection to file
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("\nSTOCH file does not exists.\n");
		return(-1);									// error: no file
	}

	// check that first line correctly declares the file as a TIME file
	if(fgets(line, sizeof line, file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field) == 1
		&& strncmp(field, "STOCH", 5) == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(longfield, problemname) == 0){
		;
	} else {
		printf("\nWrong STOCH file header.\n");
		return(-1);									// error: wrong file header
	}

	// read STOCH file type: only SCEN || SCENARIOS supported
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& (strcmp(longfield, "SCEN") == 0 || strcmp(longfield, "SCENARIOS") == 0)){
		;
	} else {
		printf("\nUnknown STOCH type declaration: %s. Only SCEN || SCENARIOS supported.\n", longfield);
		return(-1);									// error: unknown period declaration
	}

	// read line-by-line looking for " SC " keys -> number of scenarios
	*scenarionamelen = 0;
	*numscenarios = 0;
	while(fgets(line, sizeof line, file) !=  NULL){
		if(strncmp(line, " SC ", 4) == 0
			&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1){
			*scenarionamelen = (int) fmax(*scenarionamelen, strlen(longfield));
			(*numscenarios)++;
		}
	}

	// close connection with file
	fclose(file);

	// Return success indicator
	return(0);
}

/* Function to read scenario metadata */
int get_scenario_metadata(const char *problemname, const char *filename,
	const struct string_buffer *tstages,
	struct string_buffer *scen, struct string_buffer *scen_parent,
	int* branch_tstage, double *scen_prob, int *scen_start, int *scen_end)
{
	char line[SMPS_MAXLINELEN+1];
	char field[SMPS_MAXFIELDLEN+1], longfield[SMPS_MAXLINELEN+1];
	int i, j, pos;

	// open connection to file
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("\nSTOCH file does not exists.\n");
		return(-1);									// error: no file
	}

	// check that first line correctly declares the file as a TIME file
	if(fgets(line, sizeof line, file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field) == 1
		&& strcmp(field, "STOCH") == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(longfield, problemname) == 0){
		;
	} else {
		printf("\nWrong STOCH file header.\n");
		return(-1);									// error: wrong file header
	}

	// read STOCH file type: only SCEN || SCENARIOS supported
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& (strcmp(longfield, "SCEN") == 0 || strcmp(longfield, "SCENARIOS") == 0)){
		;
	} else {
		printf("\nUnknown STOCH type declaration: %s. Only SCEN || SCENARIOS supported.\n", longfield);
		return(-1);									// error: unknown STOCH declaration
	}
	
	// read line-by-line looking for " SC " keys
	i = 0;
	j = 2;
	int linelen;
	while(fgets(line, sizeof line, file) !=  NULL){
		if(strncmp(line, " SC ", 4) == 0){
			if( i >= (*scen).num_elems ){
				printf("\nSTOCH file contains at least %d scenarios, expecting only %d scenarios.\n",
					i+1, (*scen).num_elems);
				return(-1);								// error: more scenarios than expected
			}
			// get line length
			linelen = strlen(line);
			// read scenario name
			pos = SMPS_FIRSTNAMECOL;
			sscanf(&line[pos], "%" SMPS_MAXLINELENSTR "s", longfield);
			add_spaces(longfield, (*scen).atom_len+1, (*scen).atom_len,
				(*scen).buffer + ((*scen).atom_len+1)*i);
			pos = next_nonempty_pos(line, linelen, pos + strlen(longfield) - 1);
			// read scenario parent
			sscanf(&line[pos], "%" SMPS_MAXLINELENSTR "s", longfield);
			add_spaces(longfield, (*scen_parent).atom_len+1, (*scen_parent).atom_len,
				(*scen_parent).buffer + ((*scen_parent).atom_len+1)*i);
			pos = next_nonempty_pos(line, linelen, pos + strlen(longfield) - 1);
			// read scenario probability
			sscanf(&line[pos], "%" SMPS_MAXLINELENSTR "s", longfield);
			sscanf(longfield, "%lf", scen_prob + i);
			pos = next_nonempty_pos(line, linelen, pos + strlen(longfield) - 1);
			// read branching time stage
			sscanf(&line[pos], "%s", longfield);
			*(branch_tstage + i) = get_string_index(tstages, &longfield[0]);
			if( *(branch_tstage + i) < 0 ){
				printf("\nBranching period of scenario '%s' ('%s') not declared in TIME file.\n",
					(*scen).buffer + ((*scen).atom_len+1)*i, longfield);
				printf("%s\n", line);
				return(-1);							// error: unknown STOCH declaration
			}
			// start and end position within STOCH file
			*(scen_start + i) = j;
			if( i >= 1 ){
				*(scen_end + (i - 1)) = j - 1;
			}
			i++;
		}
		j++;
	}
	
	// checking I have read data for all scenarios
	if(i < (*scen).num_elems){
		printf("\nSTOCH file contains %d scenarios, expecting %d scenarios.\n",
			i, (*scen).num_elems);
		return(-1);									// error: more scenarios than expected
	}

	// end position of the last scenario
	*(scen_end + (i - 1)) = j - 2;

	// close connection with file
	fclose(file);

	// Return success indicator
	return(0);
}

/* Function to read scenario differences from the SCENARIO file */
int read_scen_differences(const char *problemname, const char *filename, const char *scenarioname,
	const struct string_buffer *scen, const struct string_buffer *scen_parent, const int *scen_start,
	const struct string_buffer *rows, const struct string_buffer *cols,
	struct scen_differences *scen_diff_values){

	int scenarioix;
	char rootname[(*scen_parent).atom_len+1];

	// Get index of scenarioname
	scenarioix = get_string_index(scen, scenarioname);
	if(scenarioix < 0){
		printf("\nScenario %s not found within scenarios buffer.\n", scenarioname);
		return(-1);
	}
	add_spaces(SMPS_ROOTNODENAME, (*scen_parent).atom_len+1, (*scen_parent).atom_len, rootname);
	if(strcmp((*scen_parent).buffer+((*scen_parent).atom_len+1)*scenarioix, rootname) != 0){
		printf("\nScenario %s has parent %s. Nested declaration of scenarios not supported.\n",
			scenarioname, (*scen_parent).buffer+((*scen_parent).atom_len+1)*scenarioix);
		return(-1);
	}

	char line[SMPS_MAXLINELEN+1];
	char code[SMPS_CODELEN+1], field1[SMPS_MAXFIELDLEN+1],
		field2[SMPS_MAXFIELDLEN+1], longfield[SMPS_MAXLINELEN+1];
	int lncnt = 0;

	// open connection to file
	FILE *file = fopen(filename, "r");
	if (file == NULL){
		printf("\nSTOCH file does not exists.\n");
		return(-1);									// error: no file
	}

	// check that first line correctly declares the file as a STOCH file
	if(fgets(line, sizeof line, file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field1) == 1
		&& strcmp(field1, "STOCH") == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(longfield, problemname) == 0){
		lncnt++;
	} else {
		printf("\nWrong STOCH file header.\n");
		return(-1);									// error: wrong file header
	}

	// read STOCH file type: only SCEN || SCENARIOS supported
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& (strcmp(longfield, "SCEN") == 0 || strcmp(longfield, "SCENARIOS") == 0)){
		lncnt++;
	} else {
		printf("\nUnknown STOCH type declaration: %s. Only SCEN || SCENARIOS supported.\n", longfield);
		return(-1);									// error: unknown STOCH declaration
	}

	// skip lines until we get to the incumbent scenario
	if( flineskip(file, *(scen_start + scenarioix)-lncnt, SMPS_MAXLINELEN+1) == 0 ){
		lncnt = *(scen_start + scenarioix);
	} else {
		printf("\nUnable to skip %d lines in STOCH file.\n", *(scen_start + scenarioix)-lncnt);
		return(-1);
	}
	
	// check scenario header coincides with requested scenario
	if(fgets(line,sizeof line,file) !=  NULL
		&& strncmp(line, " SC ", 4) == 0
		&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		// 2017-03-18: corrected for scenarios with name of varying length
		&& strcmp(scenarioname, add_spaces_in_place(longfield, (*scen).atom_len + 1, (*scen).atom_len)) == 0){
		lncnt++;
	} else {
		printf("\nLine %d of STOCH file is not coherent with a header for scenario '%s'.",
			*(scen_start + scenarioix), scenarioname);
		printf("\nlongfield: '%s'", longfield);
		return(-1);
	}

	// declaring and allocating buffers to temporally store names of modifications
	struct string_buffer rhschg_con, bndchg_var, coeffchg_var, coeffchg_con;
	int numrhschg, numbndchg, numcoeffchg;
	rhschg_con.num_elems = (*scen_diff_values).numchgrhs;
	rhschg_con.atom_len = (*rows).atom_len;
	rhschg_con.buffer = (char *) malloc(sizeof(char)*(rhschg_con.atom_len+1)*rhschg_con.num_elems);
	bndchg_var.num_elems = (*scen_diff_values).numchgbnd;
	bndchg_var.atom_len = (*cols).atom_len;
	bndchg_var.buffer = (char *) malloc(sizeof(char)*(bndchg_var.atom_len+1)*bndchg_var.num_elems);
	coeffchg_var.num_elems = (*scen_diff_values).numchgcoeff;
	coeffchg_var.atom_len = (*cols).atom_len;
	coeffchg_var.buffer = (char *) malloc(sizeof(char)*(coeffchg_var.atom_len+1)*coeffchg_var.num_elems);
	coeffchg_con.num_elems = (*scen_diff_values).numchgcoeff;
	coeffchg_con.atom_len = (*rows).atom_len;
	coeffchg_con.buffer = (char *) malloc(sizeof(char)*(coeffchg_con.atom_len+1)*coeffchg_con.num_elems);
	numrhschg = 0;
	numbndchg = 0;
	numcoeffchg = 0;
	// read scenario modifications line-by-line
	while(fgets(line,sizeof line,file) !=  NULL
		&& strncmp(line, SMPS_SCENARIO_CODE, strlen(SMPS_SCENARIO_CODE)) != 0
		&& strncmp(line, SMPS_ENDDATA, strlen(line)-1) != 0									// strlen - 1 to remove line jump
		&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field1) == 1){
		lncnt++;
		// parse line according to first field
		if( strcmp(field1, SMPS_RHSNAME) == 0 || strcmp(field1, SMPS_RHSNAME_ALT) == 0){
			// parsing an rhs change
			if( sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2) == 1
				&& insert_in_buffer(&rhschg_con, numrhschg, field2) == 0
				&& sscanf(&line[SMPS_FIRSTNUMCOL], "%lf",
					((*scen_diff_values).chgrhs_scenval + numrhschg)) == 1){
				numrhschg++;
			} else {
				printf("\nError while reading an RHS change for scenario %s from STOCH file.\n", scenarioname);
				return(-1);
			}
		} else if( strcmp(field1, SMPS_BNDNAME) == 0 || strcmp(field1, SMPS_BNDNAME_ALT) == 0 ) {
			// parsing a bound change
			if(sscanf(&line[SMPS_CODECOL], "%" SMPS_CODELENSTR "s", code) == 1
				&& (strncmp(code, "U", 1) == 0 || strncmp(code, "L", 1) == 0)
				&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2) == 1
				&& insert_in_buffer(&bndchg_var, numbndchg, field2) == 0
				&& sscanf(&line[SMPS_FIRSTNUMCOL], "%lf",
					((*scen_diff_values).chgbnd_scenval + numbndchg)) == 1){
				*((*scen_diff_values).chgbnd_type + numbndchg) = code[0];
				numbndchg++;
			} else {
				printf("\nError while reading a BND change for scenario %s from STOCH file.\n", scenarioname);
				return(-1);
			}
		} else {
			// parsing a coefficient change
			if(insert_in_buffer(&coeffchg_var, numcoeffchg, field1) == 0
				&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2) == 1
				&& insert_in_buffer(&coeffchg_con, numcoeffchg, field2) == 0
				&& sscanf(&line[SMPS_FIRSTNUMCOL], "%lf",
					((*scen_diff_values).chgcoeff_scenval + numcoeffchg)) == 1){
				numcoeffchg++;
			} else {
				printf("\nError while reading a COEFF change for scenario %s from STOCH file.\n", scenarioname);
				return(-1);
			}
		}
	}
	
	// close connection with file
	fclose(file);

	// Shrink buffers
	rhschg_con.num_elems = numrhschg;
	rhschg_con.buffer = (char *) realloc(rhschg_con.buffer,
		sizeof(char)*(rhschg_con.atom_len+1)*rhschg_con.num_elems);
	bndchg_var.num_elems = numbndchg;
	bndchg_var.buffer = (char *) realloc(bndchg_var.buffer,
		sizeof(char)*(bndchg_var.atom_len+1)*bndchg_var.num_elems);
	coeffchg_var.num_elems = numcoeffchg;
	coeffchg_var.buffer = (char *) realloc(coeffchg_var.buffer,
		sizeof(char)*(coeffchg_var.atom_len+1)*coeffchg_var.num_elems);
	coeffchg_con.num_elems = numcoeffchg;
	coeffchg_con.buffer = (char *) realloc(coeffchg_con.buffer,
		sizeof(char)*(coeffchg_con.atom_len+1)*coeffchg_con.num_elems);

	// Update number of changes
	(*scen_diff_values).numchgrhs = numrhschg;
	(*scen_diff_values).numchgbnd = numbndchg;
	(*scen_diff_values).numchgcoeff = numcoeffchg;

	// Parse modifications to column/row index format
	int i;
	int *COREsort;
	// parsing rows
	COREsort = (int *) malloc((*rows).num_elems*sizeof(int));
	for(i = 0; i < (*rows).num_elems; i++) *(COREsort + i) = i;
	sort_r(COREsort, (*rows).num_elems, sizeof(int), &compare_atoms_buffer, (void*) rows);
	match_to_presorted_buffer(rows, COREsort, &rhschg_con, (*scen_diff_values).chgrhs_row);			// rhs changes
	match_to_presorted_buffer(rows, COREsort, &coeffchg_con, (*scen_diff_values).chgcoeff_row);		// coefficient changes
	free(COREsort);
	// parsing columns
	COREsort = (int *) malloc((*cols).num_elems*sizeof(int));
	for(i = 0; i < (*cols).num_elems; i++) *(COREsort + i) = i;
	sort_r(COREsort, (*cols).num_elems, sizeof(int), &compare_atoms_buffer, (void*) cols);
	match_to_presorted_buffer(cols, COREsort, &bndchg_var, (*scen_diff_values).chgbnd_col);			// bnd changes
	match_to_presorted_buffer(cols, COREsort, &coeffchg_var, (*scen_diff_values).chgcoeff_col);		// coefficient changes
	free(COREsort);

	// Check resulting indexes -> all row and column names must have been found to continue
	int fail = 0;
	int fprev = 0;
	// rows
	for(i = 0; i < rhschg_con.num_elems; i++){
		if( *((*scen_diff_values).chgrhs_row + i) < 0 ){
			if( fail == fprev){
				printf("\nThe following row names from the STOCH file (scenario %s, RHS section),"
					"\nwere not found among CORE row names:", scenarioname);
				fail++;
			}
			printf("\n%s", rhschg_con.buffer + (rhschg_con.atom_len + 1)*i);
		}
	}
	fprev = fail;
	for(i = 0; i < coeffchg_con.num_elems; i++){
		if( *((*scen_diff_values).chgcoeff_row + i) < 0 ){
			if( fail == fprev){
				printf("\nThe following row names from the STOCH file (scenario %s, COEFF section),"
					"\nwere not found among CORE row names:", scenarioname);
				fail++;
			}
			printf("\n%s", coeffchg_con.buffer + (coeffchg_con.atom_len + 1)*i);
		}
	}
	// columns
	fprev = fail;
	for(i = 0; i < bndchg_var.num_elems; i++){
		if( *((*scen_diff_values).chgbnd_col + i) < 0 ){
			if( fail == fprev){
				printf("\nThe following column names from the STOCH file (scenario %s, BND section),"
					"\nwere not found among CORE column names:", scenarioname);
				fail++;
			}
			printf("\n%s", bndchg_var.buffer + (bndchg_var.atom_len + 1)*i);
		}
	}
	fprev = fail;
	for(i = 0; i < coeffchg_var.num_elems; i++){
		if( *((*scen_diff_values).chgcoeff_col + i) < 0 ){
			if( fail == fprev){
				printf("\nThe following column names from the STOCH file (scenario %s, COEFF section),"
					"\nwere not found among CORE column names:", scenarioname);
				fail++;
			}
			printf("\n%s", coeffchg_var.buffer + (coeffchg_var.atom_len + 1)*i);
		}
	}
	// if any row/column name was not found, return error
	if( fail > 0 ){
		printf("\nError while parsing row/column names for scenario %s from STOCH file.\n", scenarioname);
		return(-1);
	}

	// Free allocated buffers
	free(rhschg_con.buffer);
	free(bndchg_var.buffer);
	free(coeffchg_var.buffer);
	free(coeffchg_con.buffer);

	// Return success indicator
	return(0);
}

/* Function to get scenario bounds on variables from STOCH file */
int get_scenario_col_bounds(const char *problemname, const char *filename, const char *scenarioname,
	const struct string_buffer *scen, const struct string_buffer *scen_parent,
	const int *scen_start, const int *scen_end,
	const struct string_buffer *subcols, double *lowerbnd, double *upperbnd)
{
	// Get index of scenarioname
	int scenarioix = get_string_index(scen, scenarioname);
	char rootname[(*scen_parent).atom_len+1];
	if(scenarioix < 0){
		printf("\nScenario %s not found within scenarios buffer.\n", scenarioname);
		return(-1);
	}
	add_spaces(SMPS_ROOTNODENAME, (*scen_parent).atom_len+1, (*scen_parent).atom_len, rootname);
	if(strcmp((*scen_parent).buffer+((*scen_parent).atom_len+1)*scenarioix, rootname) != 0){
		printf("\nScenario %s has parent %s. Nested declaration of scenarios not supported.\n",
			scenarioname, (*scen_parent).buffer+((*scen_parent).atom_len+1)*scenarioix);
		return(-1);
	}
	
	// open connection to file
	FILE *file = fopen(filename, "r");
	if (file == NULL){
		printf("\nSTOCH file does not exists.\n");
		return(-1);									// error: no file
	}
	
	// declare variables used for reading
	char line[SMPS_MAXLINELEN+1];
	char code[SMPS_CODELEN+1], field1[SMPS_MAXFIELDLEN+1],
		field2[SMPS_MAXFIELDLEN+1], longfield[SMPS_MAXLINELEN+1];
	double newbnd;
	int lncnt = 0;
	
	// check that first line correctly declares the file as a STOCH file
	if(fgets(line, sizeof line, file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXFIELDLENSTR "s", field1) == 1
		&& strcmp(field1, "STOCH") == 0
		&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(longfield, problemname) == 0){
		lncnt++;
	} else {
		printf("\nWrong STOCH file header.\n");
		return(-1);									// error: wrong file header
	}
	
	// read STOCH file type: only SCEN || SCENARIOS supported
	if(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& (strcmp(longfield, "SCEN") == 0 || strcmp(longfield, "SCENARIOS") == 0)){
		lncnt++;
	} else {
		printf("\nUnknown STOCH type declaration: %s. Only SCEN || SCENARIOS supported.\n", longfield);
		return(-1);									// error: unknown STOCH declaration
	}
	
	// skip lines until we get to the incumbent scenario
	if( flineskip(file, *(scen_start + scenarioix)-lncnt, SMPS_MAXLINELEN+1) == 0 ){
		lncnt = *(scen_start + scenarioix);
	} else {
		printf("\nUnable to skip %d lines in STOCH file.\n", *(scen_start + scenarioix)-lncnt);
		return(-1);
	}
	
	// check scenario header coincides with requested scenario
	if(fgets(line,sizeof line,file) !=  NULL
		&& strncmp(line, " SC ", 4) == 0
		&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXLINELENSTR "s", longfield) == 1
		&& strcmp(scenarioname, add_spaces_in_place(longfield, (*scen).atom_len + 1, (*scen).atom_len)) == 0){
		lncnt++;
	} else {
		printf("\nLine %d of STOCH file is not coherent with a header for scenario '%s'.",
			*(scen_start + scenarioix), scenarioname);
		printf("\nlongfield: '%s'", longfield);
		return(-1);
	}
	
	// declaring and allocating buffers to temporally store bound modifications
	int numlbndchg = 0, numubndchg = 0;
	struct string_buffer lbndchg_var, ubndchg_var;
	lbndchg_var.num_elems = (*(scen_end + scenarioix)) - (*(scen_start + scenarioix));
	lbndchg_var.atom_len = (*subcols).atom_len;
	lbndchg_var.buffer = (char*) malloc(sizeof(char)*(lbndchg_var.atom_len+1)*lbndchg_var.num_elems);
	ubndchg_var.num_elems = (*(scen_end + scenarioix)) - (*(scen_start + scenarioix));
	ubndchg_var.atom_len = (*subcols).atom_len;
	ubndchg_var.buffer = (char*) malloc(sizeof(char)*(ubndchg_var.atom_len+1)*ubndchg_var.num_elems);
	double *lbndchg_value = (double*) malloc(sizeof(double)*lbndchg_var.num_elems);
	double *ubndchg_value = (double*) malloc(sizeof(double)*ubndchg_var.num_elems);
	
	// read scenario modifications line-by-line
	int readerror = 0;
	while(fgets(line,sizeof line,file) !=  NULL
		&& strncmp(line, SMPS_SCENARIO_CODE, strlen(SMPS_SCENARIO_CODE)) != 0
		&& strncmp(line, SMPS_ENDDATA, strlen(line)-1) != 0									// strlen - 1 to remove line jump
		&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field1) == 1){
		lncnt++;
		// parse only bound change lines (ignore the rest)
		if( strcmp(field1, SMPS_BNDNAME) == 0 || strcmp(field1, SMPS_BNDNAME_ALT) == 0 ) {
			// read the data into 
			if(sscanf(&line[SMPS_CODECOL], "%" SMPS_CODELENSTR "s", code) == 1
				&& sscanf(&line[SMPS_SECONDNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field2) == 1
				&& sscanf(&line[SMPS_FIRSTNUMCOL], "%lf", &newbnd) ){
				// lower or upper bound?
				if( code[0] == 'L' ){
					readerror = insert_in_buffer(&lbndchg_var, numlbndchg, field2);
					*(lbndchg_value + numlbndchg) = newbnd;
					numlbndchg++;
				} else if( code[0] == 'U' ) {
					readerror = insert_in_buffer(&ubndchg_var, numubndchg, field2);
					*(ubndchg_value + numubndchg) = newbnd;
					numubndchg++;
				}
				if(readerror != 0){
					printf("\nError while storing a variable name for a BND change for scenario %s from STOCH file.\n", scenarioname);
					return(-1);
				}
			} else {
				printf("\nInconsistent fields while reading a BND change for scenario %s from STOCH file.\n", scenarioname);
				return(-1);
			}
		}
	}
	
	// close connection with file
	fclose(file);
	
	// Sort interest variables alphabetically
	int *SUBsort;
	if( numlbndchg > 0 || numubndchg > 0 ){
		SUBsort = (int*) malloc(sizeof(int) * (*subcols).num_elems);
		int i;
		for(i = 0; i < (*subcols).num_elems; i++) *(SUBsort + i) = i;
		sort_r(SUBsort, (*subcols).num_elems, sizeof(int),
			&compare_atoms_buffer, (void*) subcols);
	}
	
	
	// Parse lower bound changes
	if(numlbndchg > 0){
		
		// shrink buffers
		lbndchg_var.num_elems = numlbndchg;
		lbndchg_var.buffer = (char*) realloc(lbndchg_var.buffer,
			sizeof(char)*(lbndchg_var.atom_len+1)*lbndchg_var.num_elems);
		lbndchg_value = (double*) realloc(lbndchg_value, sizeof(double)*lbndchg_var.num_elems);
		
		// match read column names to the columns of interests (subcols)
		int *lbndchg_col = (int*) malloc(sizeof(int)*lbndchg_var.num_elems);
		match_to_presorted_buffer(subcols, SUBsort, &lbndchg_var, lbndchg_col);
		
		// update bounds using the read information
		int i;
		for(i = 0; i < lbndchg_var.num_elems; i++){
			if( *(lbndchg_col + i) >= 0
				&& *(lbndchg_value + i) > *(lowerbnd + *(lbndchg_col + i)) ){
				*(lowerbnd + *(lbndchg_col + i)) = *(lbndchg_value + i);
			}
		}
		
		// free space allocated within the scope
		free(lbndchg_col);
		
	}
	
	// Parse upper bound changes
	if(numubndchg > 0){
		
		// shrink buffers
		ubndchg_var.num_elems = numubndchg;
		ubndchg_var.buffer = (char*) realloc(ubndchg_var.buffer,
			sizeof(char)*(ubndchg_var.atom_len+1)*ubndchg_var.num_elems);
		ubndchg_value = (double*) realloc(ubndchg_value, sizeof(double)*ubndchg_var.num_elems);
		
		// match read column names to the columns of interests (subcols)
		int *ubndchg_col = (int*) malloc(sizeof(int)*ubndchg_var.num_elems);
		match_to_presorted_buffer(subcols, SUBsort, &ubndchg_var, ubndchg_col);
		
		// update bounds using the read informatio
		int i;
		for(i = 0; i < ubndchg_var.num_elems; i++){
			if( *(ubndchg_col + i) >= 0
				&& *(ubndchg_value + i) < *(upperbnd + *(ubndchg_col + i)) ){
				*(upperbnd + *(ubndchg_col + i)) = *(ubndchg_value + i);
			}
		}
		
		// free space allocated within the scope
		free(ubndchg_col);
		
	}
	
	// Free allocated buffers
	free(lbndchg_var.buffer);
	free(ubndchg_var.buffer);
	free(lbndchg_value);
	free(ubndchg_value);
	if( numlbndchg > 0 || numubndchg > 0 ){
		free(SUBsort);
	}
	
	// Return success indicator
	return(0);
}

// Read delayed rows from CORE file
int read_delayed_rows(const char *filename, const struct string_buffer *rows,
	int *num_delayed_rows, int *delayed_rows_start)
{
	// variables for holding lines as they are read
	char line[SMPS_MAXLINELEN+1];
	char longfield[SMPS_MAXLINELEN+1];
	char field[SMPS_MAXFIELDLEN+1];

	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("\nCORE file does not exist.\n");
		return(-1);									// error: file does not exist
	}

	// Skip lines until finding DELAYEDROWS or LAZYCONS (or EOF)
	int found = 0;
	while(fgets(line,sizeof line,file) !=  NULL
		&& sscanf(line, "%" SMPS_MAXLINELENSTR "s", longfield) == 1 ){
		if( strcmp(longfield, "DELAYEDROWS") == 0 || strcmp(longfield, "LAZYCONS") == 0 ){
			found = 1;
			break;
		}
	}

	// Check wether there are delayed rows
	if( found == 0 ){
		// there are no delayed rows, return inmediately
		*num_delayed_rows = 0;
		return(0);
	}

	// Collect delayed rows
	int i = 0;
	while(fgets(line,sizeof line,file) !=  NULL
		&& line[0] == ' '
		&& sscanf(&line[SMPS_FIRSTNAMECOL], "%" SMPS_MAXFIELDLENSTR "s", field) == 1 ){
		if( i == 0 ){		// get start index
			*delayed_rows_start = get_string_index(rows, field);
			if( *delayed_rows_start < 0 ){
				printf("\nRow %s in CORE file not found among rows read by Xpress.\n", field);
				return(-1);
			}
			i++;
		}else{
			if( strcmp(add_spaces_in_place(field, SMPS_MAXFIELDLEN + 1, (*rows).atom_len),
					(*rows).buffer + ((*rows).atom_len + 1)*((*delayed_rows_start)+i)) == 0 ){
				i++;
			}else{
				printf("\nRow %s in CORE file not found in a contiguous range"
					"\nof delayed rows among rows read by Xpress.\n", field);
				return(-1);
			}
		}
	}
	*num_delayed_rows = i;

	// close connection with file
	fclose(file);

	// Return success indicator
	return(0);
}