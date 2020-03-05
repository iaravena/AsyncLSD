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
 * Functions for generating stage subproblems using Xpress C API
 */

/* SMPS Xpress algorithm header */
#include "SMPS_Xpress.h"

/* Function to get nonzero elements in a set of columns or rows */
int num_submatrix_nonzeros(const XPRSprob *scenprob, const int accessdim, const int numvecs, int *vecind)
{
	// sort provided indices in ascending order (in place)
	qsort((void*) vecind, numvecs, sizeof(int), &compare_integers);
	
	// assign functions to get nonzero elements
	int (*getelemfromXpress)(XPRSprob, int*, int*, double*, const int, int*, const int, const int);
	if(accessdim == 1){				// vecind corresponds to row numbers
		getelemfromXpress = &XPRSgetrows;
	} else if(accessdim == 2){		// vecind correspond to col numbers
		getelemfromXpress = &XPRSgetcols;
	} else {
		printf("Unrecorgnized dimension value %d (must be either 1 for rows or 2 for columns).", accessdim);
		return(-1);
	}
	
	// get number of nonzero elements from Xpress
	int i = 0, j, retcoeffs, cumcoeffs = 0;
	while( i < numvecs){
		// detect whether there is a contiguous range ahead of i
		j = 1;
		while( (*(vecind + i + j)) - (*(vecind + i)) == j ) j++;
		// get coefficients from Xpress
		getelemfromXpress(*scenprob, NULL, NULL, NULL, 0, &retcoeffs, *(vecind + i), *(vecind + i + j - 1));
		// go to the next contiguous range
		cumcoeffs = cumcoeffs + retcoeffs;
		i = i + j;
	}
	
	// retrieve number of coefficients from Xpress
	return(cumcoeffs);
}

/* Function to get nonzero elements in a set of columns or rows */
int get_nonzero_elems(const XPRSprob *scenprob,
	int *mstart, int *mxind, double *dmatval, 
	const size_t maxretcoeffs, int *numcoeffs,
	const int accessdim, const int numind, int *vecind)
{
	// assign functions to get nonzero elements
	int (*getelemfromXpress)(XPRSprob, int*, int*, double*, const int, int*, const int, const int);
	if(accessdim == 1){				// vecind corresponds to row numbers
		getelemfromXpress = &XPRSgetrows;
	} else if(accessdim == 2){		// vecind correspond to col numbers
		getelemfromXpress = &XPRSgetcols;
	} else {
		printf("\nUnrecorgnized dimension value %d (must be either 1 for rows or 2 for columns).\n", accessdim);
		return(-1);
	}
	
	// sort provided indices in ascending order (in place)
	qsort((void*) vecind, numind, sizeof(int), &compare_integers);
	
	// retrieve coefficients from Xpress
	int i, j, k, retcoeffs = 0, cumcoeffs = 0;
	int *mstartlocal=NULL, *mxindlocal=NULL;
	double *dmatvallocal=NULL;
	if(mstart != NULL){
		mstartlocal = (int*) malloc(sizeof(int)*(numind+1));
		*mstart = 0;
	}
	i = 0;
	while( i < numind && cumcoeffs < maxretcoeffs ){
		// detect whether there is a contiguous range ahead of i
		j = 1;
		while( (*(vecind + i + j)) - (*(vecind + i)) == j ) j++;
		// get coefficients from Xpress
		if(mxind != NULL) mxindlocal = mxind + cumcoeffs;
		if(dmatval != NULL) dmatvallocal = dmatval + cumcoeffs;
		getelemfromXpress(*scenprob, mstartlocal, mxindlocal, dmatvallocal,
			maxretcoeffs - cumcoeffs, &retcoeffs, *(vecind + i) , *(vecind + i + j - 1));
		// fill mstart
		if( mstart != NULL ){
			for(k = 1; k <= j; k++)
				*(mstart + i + k) = (*(mstart + i)) + (*(mstartlocal + k));
		}
		// go to the next contiguous range
		cumcoeffs = cumcoeffs + retcoeffs;
		i = i + j;
	}
	*numcoeffs = cumcoeffs;
	if(mstart != NULL) free(mstartlocal);
	
	// sort minor (e.g. sort rows if we obtained the value by columns)
	if( mstart != NULL && mxind != NULL){
		int *auxind, blockstart, blocksize;
		auxind = (int*) malloc(sizeof(int)*(*numcoeffs));
		for(i = 0; i < *numcoeffs; i++) *(auxind + i) = i;
		for(i = 0; i < numind; i++){
			blockstart = *(mstart + i);
			blocksize = (*(mstart + i + 1)) - (*(mstart + i));
			sort_r((void*) (auxind + blockstart), blocksize, sizeof(int), &compare_integers_r, (void*) mxind);
		}
		if(dmatval != NULL){		// order mxind directly
			reorder_int_vector(*numcoeffs, mxind, auxind);
		}else{						// generate a data frame to order mxind and dmatval together
			struct data_frame auxdf;
			auxdf.numcols = 2;
			auxdf.numrows = *numcoeffs;
			auxdf.coltypes = (int*) malloc(sizeof(int)*auxdf.numcols);
			*(auxdf.coltypes) = DFR_COLTYPE_INT;
			*(auxdf.coltypes + 1) = DFR_COLTYPE_DOUBLE;
			auxdf.ptr = (void**) malloc(sizeof(void*)*auxdf.numcols);
			*(auxdf.ptr) = (void*) mxind;
			*(auxdf.ptr + 1) = (void*) dmatval;
			reorder_data_frame(&auxdf, auxind);
			free(auxdf.coltypes);
			free(auxdf.ptr);
		}
		free(auxind);
	}
	
	// return success indicator
	return(0);
}

/* Function to transpose matrix arrays returned by Xpress */
int transpose_sparse_matrix(
	const int in_numind, const int *in_vecind, const int* in_mstart,				// original major order
	const int out_numind, const int *out_vecind, int *out_mstart,					// original minor order
	const int numnonzeros, const int *in_mxind, int *out_mxind,
		const double *in_dmatval, double *out_dmatval)								// secondary indexes and values
{
	// Iterators
	int i, j, k;
	
	// Allocate vector with positions within each block
	int *pos = (int*) malloc(sizeof(int)*in_numind);
	for(j = 0; j < in_numind; j++) *(pos + j) = *(in_mstart + j);
	
	// Perform the transposition
	// (take advantage of the fact that the input is already sorted)
	*out_mstart = 0;
	for(i = 0; i < out_numind; i++){
		k = *(out_mstart + i);
		for(j = 0; j < in_numind; j++){
			if( *(pos + j) < *(in_mstart + j + 1)
				&& *(out_vecind + i) == *(in_mxind + *(pos + j)) ){
				*(out_mxind + k) = *(in_vecind + j);
				*(out_dmatval + k) = *(in_dmatval + *(pos + j));
				*(pos + j) += 1;
				k++;
			}else{
				continue;
			}
		}
		*(out_mstart + i + 1) = k;
	}
	
	// Check that all values have been included in the transposition
	j = 0;
	while( j < in_numind && *(pos + j) == *(in_mstart + j + 1) ) j++;
	if( *(out_mstart + out_numind) != numnonzeros || j != in_numind ){
		printf("\nCollected elements are not consisten with input data (%d:%d,%d:%d).\n",
			*(out_mstart + out_numind + 1), numnonzeros, j, in_numind);
		return(-1);
	}
	
	// Free allocated memory
	free(pos);
	
	// return success indicator
	return(0);
}

/* Function to get row and column indices associated with a time stage */
// (rows with variables from other periods are not included in the output)
int select_cols_and_rows(const XPRSprob *sourceprob,
	const struct string_buffer *rows, const struct string_buffer *cols,
	const int num_ignore_rows, const int ignore_rows_start,
	const struct string_buffer *timestages, const int *rowtimestage, const int *coltimestage,
	const int numcoltss, int *coltss, const int numrowtss, int *rowtss,
	int *numinccols, int *inccols, int *numincrows, int *incrows, int *numnonzerosub)
{
	// iterators
	int i, j;
	
	// get columns of the indicated periods
	if(numinccols != NULL){
		int *clind = malloc(sizeof(int)*(*cols).num_elems);
		for(i = 0; i < (*cols).num_elems; i++) *(clind + i) = i;
		get_items_by_property_value((*cols).num_elems, clind, coltimestage,
			numcoltss, coltss, (*cols).num_elems, numinccols, inccols);
		free(clind);
	}
	
	// get rows (not ignored) of the indicated periods
	if(numincrows != NULL){
		int *rwind = (int*) malloc(sizeof(int)*((*rows).num_elems-num_ignore_rows));
		if( num_ignore_rows == 0 ){
			for(i = 0; i < (*rows).num_elems; i++){
				*(rwind + i) = i;
			}
			get_items_by_property_value((*rows).num_elems, rwind, rowtimestage,
				numrowtss, rowtss, (*rows).num_elems, numincrows, incrows);
		}else if( num_ignore_rows > 0 && num_ignore_rows <= (*rows).num_elems
			&& ignore_rows_start >= 0 && ignore_rows_start < (*rows).num_elems ){
			int *rtstage = (int*) malloc(sizeof(int)*((*rows).num_elems-num_ignore_rows));
			for(i = 0; i < ignore_rows_start; i++){
				*(rwind + i) = i;
				*(rtstage + i) = *(rowtimestage + i);
			}
			for(i=ignore_rows_start + num_ignore_rows; i < (*rows).num_elems; i++){
				*(rwind + i - num_ignore_rows) = i;
				*(rtstage + i - num_ignore_rows) = *(rowtimestage + i);
			}
			get_items_by_property_value((*rows).num_elems-num_ignore_rows, rwind, rtstage,
				numrowtss, rowtss, (*rows).num_elems, numincrows, incrows);
			free(rtstage);
		}else{
			printf("\nInconsistent values for ignoring rows:"
				"\nNum. rows: %d\tNum. ignore: %d\tStart ignore: %d\n",
				(*rows).num_elems, num_ignore_rows, ignore_rows_start);
			return(-1);
		}
		free(rwind);
	}
	
	// Processing rows: remove rows with non zero elements from periods outside the specified column periods
	if(numinccols != NULL && numincrows != NULL){
		
		// get number of nonzero elements in incrows and allocate buffers for get nonzero coordinates
		int expected_coeffs = num_submatrix_nonzeros(sourceprob, 1, *numincrows, incrows);
		int *mstart = (int*) malloc(sizeof(int)*((*numincrows) + 1));
		int *clind = (int*) malloc(sizeof(int)*expected_coeffs);
		int *taburows = (int*) malloc(sizeof(int)*(*numincrows));
		
		// get position of nonzero elements in incrows
		int numretcoeffs;
		if( get_nonzero_elems(sourceprob, mstart, clind, NULL,
			expected_coeffs, &numretcoeffs, 1, *numincrows, incrows) != 0){
			printf("\nError while getting non-zero elements for row cleaning.\n");
			return(-1);
		}else if( numretcoeffs != expected_coeffs){
			printf("\nNumber of elements returned by Xpress (%d) differs from the expected quantity (%d).\n",
				expected_coeffs, numretcoeffs);
			return(-1);
		}
		
		// detect rows with nonzero elements outside the column periods
		int blockstart, blocklength;
		for(i = 0; i < *numincrows; i++){
			// get same row block
			blockstart = *(mstart + i);
			blocklength = (*(mstart + i + 1)) - (*(mstart + i));
			// determine whether all columns in the row are included
			if(blocklength > 0){
				*(taburows + i) = 1 - integer_set_inclusion(blocklength,
					clind + blockstart, *numinccols, inccols);
			}else{
				*(taburows + i) = 1;
			}
		}
		
		// remove rows with nonzero elements outside the column periods from incrows, in place
		j = 0;
		for(i = 0; i < *numincrows; i++){
			if(*(taburows + i) == 0){
				if( i > j ){
					*(incrows + j) = *(incrows + i);
				}
				j++;
			}
		}
		*numincrows = j;
		
		// return number of coefficients
		if(numnonzerosub != NULL) *numnonzerosub = expected_coeffs;
		
		// free unused space
		free(mstart);
		free(clind);
		free(taburows);
	}
	
	// return success indicator
	return(0);
}

/* Function for preparing arrays to call the XPRSload_ functions with a subproblem of *sourceprob */
// extracting only MILP information from *sourceprob: inccols and incrows
int prepare_XPRSload_call(const XPRSprob *sourceprob,
	const int numinccols, const int *inccols,					// inccols must be sorted in increasing order
	const int numincrows, int *incrows,							// incrows must be sorted in increasing order
	const int numnonzerosub, const int ngentsub,				// bounds used to malloc the output buffers (this function is not safe)
	char *qrtype, double *rhs, double *obj,
	int *numnonzeros, int *mstart, int *mnel, int *mrwind, double *dmatval,
	double *dlb, double *dub, int *ngents, char *qgtype, int *mgcols)
{
	// iterators
	int i, j, k;
	
	// row types
	{
		char *qrtypeall;
		qrtypeall = (char *) malloc(sizeof(char)*( (*(incrows + numincrows - 1)) - (*incrows) + 1 ));
		XPRSgetrowtype(*sourceprob, qrtypeall, *incrows, *(incrows + numincrows - 1));
		for(i = 0; i < numincrows; i++) *(qrtype + i) = *(qrtypeall + (*(incrows + i)) - (*incrows));
		free(qrtypeall);
	}
	
	// rhs
	{
		double *rhsall;
		rhsall = (double *) malloc(sizeof(double)*( (*(incrows + numincrows - 1)) - (*incrows) + 1 ));
		XPRSgetrhs(*sourceprob, rhsall, *incrows, *(incrows + numincrows - 1));
		for(i = 0; i < numincrows; i++) *(rhs + i) = *(rhsall + (*(incrows + i)) - (*incrows));
		free(rhsall);
	}
	
	// objective
	if(obj != NULL){
		double *objall;
		objall = (double *) malloc(sizeof(double)*( (*(inccols + numinccols - 1)) - (*inccols) + 1 ));
		XPRSgetobj(*sourceprob, objall, *inccols, *(inccols + numinccols - 1));
		for(i = 0; i < numinccols; i++) *(obj + i) = *(objall + (*(inccols + i)) - (*inccols));
		free(objall);
	}
	
	// coefficients
	{
		// get the coefficients row-by-row
		int *mstartlocal = (int*) malloc(sizeof(int)*(numincrows + 1));
		int *clindlocal = (int*) malloc(sizeof(int)*numnonzerosub);
		double *dmatvallocal = (double*) malloc(sizeof(double)*numnonzerosub);
		if( get_nonzero_elems(sourceprob, mstartlocal, clindlocal, dmatvallocal,
			numnonzerosub, numnonzeros, 1, numincrows, incrows) != 0){
			printf("\nError while getting non-zero elements for preparing Xpress call.\n");
			return(-1);
		}else if( *numnonzeros > numnonzerosub ){
			printf("\nNumber of elements returned by Xpress (%d) is larger than expected (%d).\n",
				*numnonzeros, numnonzerosub);
			return(-1);
		}
		clindlocal = (int*) realloc(clindlocal, sizeof(int)*(*numnonzeros));
		dmatvallocal = (double*) realloc(dmatvallocal, sizeof(double)*(*numnonzeros));
		
		// transpose coefficients to column-wise representation
		int *relincrows = (int*) malloc(sizeof(int)*numincrows);
		for(i = 0; i < numincrows; i++) *(relincrows + i) = i;
		if( transpose_sparse_matrix(
			numincrows, relincrows, mstartlocal, numinccols, inccols, mstart,
			*numnonzeros, clindlocal, mrwind, dmatvallocal, dmatval) != 0 ){
			printf("\nProblem while transposing sparse array for preparing Xpress call.\n");
			return(-1);
		}
		
		// compute mnel
		for(i = 0; i < numinccols; i++)
			*(mnel + i) = (*(mstart + i + 1)) - (*(mstart + i));
		
		// free allocated memory
		free(relincrows);
		free(mstartlocal);
		free(clindlocal);
		free(dmatvallocal);
	}
	
	// bounds
	{
		double *bndall;
		bndall = (double *) malloc(sizeof(double)*( (*(inccols + numinccols - 1)) - (*inccols) + 1 ));
		XPRSgetlb(*sourceprob, bndall, *inccols, *(inccols + numinccols - 1));
		for(i = 0; i < numinccols; i++) *(dlb + i) = *(bndall + (*(inccols + i)) - (*inccols));
		XPRSgetub(*sourceprob, bndall, *inccols, *(inccols + numinccols - 1));
		for(i = 0; i < numinccols; i++) *(dub + i) = *(bndall + (*(inccols + i)) - (*inccols));
		free(bndall);
	}
	
	// globals (integers)
	if(ngentsub > 0){
		int nsets, *mgcolsall;
		char *qgtypeall;
		qgtypeall = (char*) malloc(sizeof(char)*ngentsub);
		mgcolsall = (int*) malloc(sizeof(int)*ngentsub);
		XPRSgetglobal(*sourceprob, ngents, &nsets, qgtypeall, mgcolsall,
			NULL, NULL, NULL, NULL, NULL);
		// sort qgtypeall and mgcolsall -> are these returned in order?? Ans: they are for our instances -> do nothing for now, but be careful
		// get positions of globals within inccols, form qgtype and mgcols
		j = 0;
		k = 0;
		for(i = 0; i < *ngents; i++){
			// move j until we find the current column
			while( *(mgcolsall + i) > *(inccols + j) && j < numinccols - 1 ) j++;
			// check whether this global belongs to inccols
			if( *(mgcolsall + i) == *(inccols + j) ){
				*(qgtype + k) = *(qgtypeall + i);
				*(mgcols + k) = j;
				k++;
			}
		}
		*ngents = k;
		free(qgtypeall);
		free(mgcolsall);
	}
	
	// return success indicator
	return(0);
}

/* Function to detect presence in multiple periods of first stage variables */
int get_firststagecols_multiplicity(const XPRSprob *scenprob,
	const struct string_buffer *rows, const struct string_buffer *cols,
	const struct string_buffer *timestages, const int *rowtimestage, const int *coltimestage,
	const int firststageind, const int num_firststgcols, const int *firststgcols,
	int *firststgcols_multiplicity)
{
	// Main loop over all first stage columns
	int i, j, k, *incrows, retrows = 0, *colcoeffsts, colnumts = 0;
	incrows = (int*) calloc((*rows).num_elems, sizeof(int));
	colcoeffsts = (int*) calloc((*timestages).num_elems, sizeof(int));
	for(i = 0; i < num_firststgcols; i++){
		// get row indices from Xpress
		XPRSgetcols(*scenprob, NULL, incrows, NULL,
			(*rows).num_elems, &retrows, *(firststgcols + i), *(firststgcols + i));
		// get unique row time stages
		for(j = 0; j < retrows; j++){
			k = 0;
			while(k < colnumts && 
				( *(rowtimestage + *(incrows + j)) != *(colcoeffsts + k) || k != firststageind )) k++;
			if(k == colnumts && colnumts < (*timestages).num_elems){
				*(colcoeffsts + colnumts) = *(rowtimestage + *(incrows + j));
				colnumts++;
			} else if(colnumts > (*timestages).num_elems){
				printf("\nThere are more time stages in column %s than there time stages in the problem.\n",
					(*cols).buffer + ((*cols).atom_len + 1) * *(firststgcols + i)) ;
				return(-1);
			}
		}
		// store multiplicity
		*(firststgcols_multiplicity + i) = colnumts;
		// erase coefficients and coefficient time stages
		for(j = 0; j < retrows; j++) *(incrows + j) = 0;
		retrows = 0;
		for(j = 0; j < colnumts; j++) *(colcoeffsts + j) = 0; 
		colnumts = 0;
	}
	free(incrows);
	free(colcoeffsts);
	
	// Return success indicator
	return(0);
}

/* Function to create scenario-period subproblem */
int load_params_period_subprob(const XPRSprob *scenprob, const char *periodprobname,
	const struct string_buffer *rows, const struct string_buffer *cols,
	const int num_ignore_rows, const int ignore_rows_start,
	const struct string_buffer *timestages, const int *rowtimestage, const int *coltimestage,
	const int firststageind, const int period,
	const int num_firststgcols, const int *firststgcols, const int *firststgcols_multiplicity,
	// first stage vars in *scenprob: number of first stage variables, first stage variable indices in *scenprob, first stage variables on multiple periods
	XPRSprob *periodprob, int *num_perfirststgcols, int *perfirststgcols, int *perfirststgcol_ind)
	// first stage vars in *periodprob: number of first stage variables, first stage variable indices in *periodprob, first stage variables indices in firststgcols
{
	// Iterators
	int i, j;
	
	// Check that $period does not corresponds to first stage
	if(period == firststageind){
		printf("\nPeriod subproblem cannot be formulated for first stage.\n");
		return(-1);
	}
	
	// Get incumbent columns, rows and an upper bound on the number of non-zero coefficients
	int colstgs[] = {firststageind, period};
	int rowstgs[] = {period};
	int numinccols, *inccols, numincrows, *incrows, numcoeffub;
	inccols = (int*) malloc(sizeof(int)*(*cols).num_elems);		// columns of the indicated 'column stages'
	incrows = (int*) malloc(sizeof(int)*(*rows).num_elems);		// rows of the indicated 'rows stages' with nonzeros only for inccols
	if( select_cols_and_rows(scenprob, rows, cols, num_ignore_rows, ignore_rows_start,
		timestages, rowtimestage, coltimestage, 2, &colstgs[0], 1, &rowstgs[0],
		&numinccols, inccols, &numincrows, incrows, &numcoeffub) != 0 ){
		printf("\nProblem detected while isolating incumbent rows and columns for period %d.\n", period);
		return(-1);
	}
	inccols = (int*) realloc(inccols, sizeof(int)*numinccols);
	incrows = (int*) realloc(incrows, sizeof(int)*numincrows);
	
	// Construct call to generate subproblem
	int ngentsub=0, numcoeffs=0, *mstart, *mnel, *mrwind, ngents=0, *mgcols;
	char *qrtype, *qgtype;
	double *rhs, *obj, *dmatval, *dlb, *dub;
	// matrix data
	qrtype = (char*) malloc(sizeof(char)*numincrows);
	rhs = (double*) malloc(sizeof(double)*numincrows);
	obj = (double*) malloc(sizeof(double)*numinccols);
	mstart = (int*) malloc(sizeof(int)*(numinccols+1));			// extra slot required to transpose without mnel
	mnel = (int*) malloc(sizeof(int)*numinccols);
	mrwind = (int*) malloc(sizeof(int)*numcoeffub);
	dmatval = (double*) malloc(sizeof(double)*numcoeffub);
	dlb = (double*) malloc(sizeof(double)*numinccols);
	dub = (double*) malloc(sizeof(double)*numinccols);
	// integer data
	XPRSgetintattrib(*scenprob, XPRS_MIPENTS, &ngentsub);
	qgtype = (char*) malloc(sizeof(char)*ngentsub);
	mgcols = (int*) malloc(sizeof(int)*ngentsub);
	// prepare call
	if( prepare_XPRSload_call(scenprob,
		numinccols, inccols, numincrows, incrows,
		numcoeffub, ngentsub, qrtype, rhs, obj,
		&numcoeffs, mstart, mnel, mrwind, dmatval,
		dlb, dub, &ngents, qgtype, mgcols) != 0){
		printf("\nError while preparing call to XPRSload for generating period subproblem.\n");
		return(-1);
	}
	// resize buffers to their actual size
	mstart = (int*) realloc(mstart, sizeof(int)*numinccols);
	mrwind = (int*) realloc(mrwind, sizeof(int)*numcoeffs);
	dmatval = (double*) realloc(dmatval, sizeof(double)*numcoeffs);
	qgtype = (char*) realloc(qgtype, sizeof(char)*ngents);
	mgcols = (int*) realloc(mgcols, sizeof(int)*ngents);
	
	// Find first stage cols in period subproblem
	int firststageint = firststageind;
	int *perindcol = malloc(sizeof(int)*numinccols);
	int *inccolts = malloc(sizeof(int)*numinccols);
	for(i = 0; i < numinccols; i++) *(perindcol + i) = i;
	for(i = 0; i < numinccols; i++) *(inccolts + i) = *(coltimestage + *(inccols + i));
	get_items_by_property_value(numinccols, perindcol, inccolts, 1, &firststageint,
		num_firststgcols, num_perfirststgcols, perfirststgcols);
	j = 0;
	for(i = 0; i < *num_perfirststgcols; i++){
		while( j < num_firststgcols &&
			*(inccols + *(perfirststgcols + i)) > *(firststgcols + j) ) j++;
		if( *(inccols + *(perfirststgcols + i)) == *(firststgcols + j) ){
			*(perfirststgcol_ind + i) = j;
		} else {
			printf("\nError found while detecting first stage columns in period subproblem.\n");
			//printf("\nHey!\n%d\t%d\t%d\t%d\n",
			//	i, j, *(inccols + *(perfirststgcols + i)), *(firststgcols + j));
			//for(j = 0; j < *num_perfirststgcols; j++)
			//	printf("\n%d,%d\t", *(perfirststgcols + j), *(inccols + *(perfirststgcols + j)));
			//for(i = 0; i < numinccols; i++){
			//	printf("\n%d\t%d\t%d", *(perindcol + i), *(inccols + i), *(inccolts + i));
			//}
			return(-1);
		}
	}
	free(perindcol);
	free(inccolts);
	
	// Scale objective coefficients of first stage columns according
	// to their period multiplicity (firststgcols_multiplicity)
	j = 0;
	for(i = 0; i < *num_perfirststgcols; i++){
		if( *(obj + *(perfirststgcols + i)) != 0
			&& *(firststgcols_multiplicity + *(perfirststgcol_ind + i)) > 1){
			*(obj + *(perfirststgcols + i)) = *(obj + *(perfirststgcols + i)) /
				((double) *(firststgcols_multiplicity + *(perfirststgcol_ind + i)) );
		}
	}
	
	// Load MILP into the optimizer
	XPRSloadglobal(*periodprob, periodprobname, numinccols, numincrows,
		qrtype, rhs, NULL, obj, mstart, mnel, mrwind, dmatval, dlb, dub,
		ngents, 0, qgtype, mgcols, NULL, NULL, NULL, NULL, NULL);
	
	// Free heap
	free(inccols);
	free(incrows);
	free(qrtype);
	free(rhs);
	free(obj);
	free(mstart);
	free(mnel);
	free(mrwind);
	free(dmatval);
	free(dlb);
	free(dub);
	free(qgtype);
	free(mgcols);
	
	// Return success indicator
	return(0);
}

/* Function to create a subproblem containing only the first stage matrix */
int load_first_stage_matrix(const XPRSprob *coreprob,
	const struct string_buffer *rows, const int *rowtimestage,
	const int firststageind, const int num_firststgcols, const int *firststgcols,
	XPRSprob *firststageprob, const char *firststageprobname)
{
	// Iterators
	int i;
	
	// Rows variables
	int numincrows, *incrows;
	incrows = (int*) malloc(sizeof(int)*(*rows).num_elems);
	
	// Get incumbent rows: here we take a leap of faith, we do not delete rows marked as first period that have
	// non zeros for columns outside the first stage columns
	{
		int *rwind = malloc(sizeof(int)*(*rows).num_elems);
		for(i = 0; i < (*rows).num_elems; i++) *(rwind + i) = i;
		int rowtss[] = {firststageind};
		get_items_by_property_value((*rows).num_elems, rwind, rowtimestage,
			1, &rowtss[0], (*rows).num_elems, &numincrows, incrows);
		free(rwind);
		incrows = (int*) realloc(incrows, sizeof(int)*numincrows);
	}
	
	// Generate arrays for calling Xpress
	int numcoeffub = num_submatrix_nonzeros(coreprob, 1, numincrows, incrows), ngentsub = 0;
	int numcoeffs = 0, *mstart, *mnel, *mrwind, ngents = 0, *mgcols;
	char *qrtype, *qgtype;
	double *rhs, *obj, *dmatval, *dlb, *dub;
	// matrix data
	qrtype = (char*) malloc(sizeof(char)*numincrows);
	rhs = (double*) malloc(sizeof(double)*numincrows);
	obj = (double*) malloc(sizeof(double)*num_firststgcols);
	for(i = 0; i < num_firststgcols; i++) *(obj + i) = 0;
	mstart = (int*) malloc(sizeof(int)*(num_firststgcols+1));			// extra slot required to transpose without mnel
	mnel = (int*) malloc(sizeof(int)*num_firststgcols);
	mrwind = (int*) malloc(sizeof(int)*numcoeffub);
	dmatval = (double*) malloc(sizeof(double)*numcoeffub);
	dlb = (double*) malloc(sizeof(double)*num_firststgcols);
	dub = (double*) malloc(sizeof(double)*num_firststgcols);
	// integer data
	XPRSgetintattrib(*coreprob, XPRS_MIPENTS, &ngentsub);
	qgtype = (char*) malloc(sizeof(char)*ngentsub);
	mgcols = (int*) malloc(sizeof(int)*ngentsub);
	// prepare call
	if( prepare_XPRSload_call(coreprob,
		num_firststgcols, firststgcols, numincrows, incrows,
		numcoeffub, ngentsub, qrtype, rhs, NULL,
		&numcoeffs, mstart, mnel, mrwind, dmatval,
		dlb, dub, &ngents, qgtype, mgcols) != 0){
		printf("\nError while preparing call to XPRSload for loading first stage matrix.\n");
		return(-1);
	}
	// resize buffers to their actual size
	mstart = (int*) realloc(mstart, sizeof(int)*num_firststgcols);
	mrwind = (int*) realloc(mrwind, sizeof(int)*numcoeffs);
	dmatval = (double*) realloc(dmatval, sizeof(double)*numcoeffs);
	qgtype = (char*) realloc(qgtype, sizeof(char)*ngents);
	mgcols = (int*) realloc(mgcols, sizeof(int)*ngents);
	
	// Load problem as a MIQP
	XPRSloadqglobal(*firststageprob, firststageprobname,
		num_firststgcols, numincrows, qrtype, rhs, NULL, obj,
		mstart, mnel, mrwind, dmatval, dlb, dub, 0, NULL, NULL, NULL,
		ngents, 0, qgtype, mgcols, NULL, NULL, NULL, NULL, NULL);
	
	// Free heap
	free(incrows);
	free(qrtype);
	free(rhs);
	free(mstart);
	free(mnel);
	free(mrwind);
	free(dmatval);
	free(dlb);
	free(dub);
	free(qgtype);
	free(mgcols);
	
	// Return success indicator
	return(0);
}

/* Function to load objective for non-scenario MIQP/QP subproblem:
 * f_0^mu(-sum x) = min_{u in U} {((-sum x)./(u_UB - u_LB))^T u + .5*mu*|(u - u0)./(u_UB - u_LB)|^2} */
int load_f0_quad_objective(const int simplifybinaries, const struct string_buffer *cols,
	const int num_firststgcols, const int *firststgcols,
	const double *firststagelb, const double *firststageub,
	const double *msumx, const double mu, const double *u0,
	XPRSprob *firststageprob, double *objoffset)					// firststageprob should already have the constraint matrix
{
	// Iterators
	int i, j;
	
	// Check that problem dimensions agree with firts stage variables
	int ncolsxprs;
	XPRSgetintattrib(*firststageprob, XPRS_ORIGINALCOLS, &ncolsxprs);
	if( ncolsxprs != num_firststgcols ){
		printf("\nProvided XPRSprob (firststageprob) has a number of columns (%d) incompatible"
			"\nwith the provided number of first stage columns (%d).\n\n",
			ncolsxprs, num_firststgcols);
		return(-1);
	}
	
	// declare quadratic objective arrays
	char *coltype = (char*) malloc(sizeof(char)*ncolsxprs);
	XPRSgetcoltype(*firststageprob, coltype, 0, ncolsxprs-1);
	double *obj, *dqe;
	int nqtr, *mqc1, *mqc2;
	obj = (double*) malloc(sizeof(double)*num_firststgcols);
	if(simplifybinaries == 1){
		XPRSgetintattrib(*firststageprob, XPRS_MIPENTS, &nqtr);
		nqtr = ncolsxprs - nqtr;
	} else {
		nqtr = num_firststgcols;
	}
	mqc1 = (int*) malloc(sizeof(int)*nqtr);
	mqc2 = (int*) malloc(sizeof(int)*nqtr);
	dqe = (double*) malloc(sizeof(double)*nqtr);
	
	// compute linear and quadratic scaling factors
	double *eta = (double*) malloc(sizeof(double)*ncolsxprs);
	double *kappa = (double*) malloc(sizeof(double)*ncolsxprs);
	for(i = 0; i < ncolsxprs; i++){
		if( (*(firststageub + i)) - (*(firststagelb + i)) > DBL_EPSILON 
			&& *(coltype + i) != 'B' ){
			*(eta + i) = 1.0/((*(firststageub + i)) - (*(firststagelb + i)));
			*(kappa + i) = (*(eta + i)) * (*(eta + i));
		}else{
			*(eta + i) = 1.0;
			*(kappa + i) = 1.0;
		}
	}
	
	// compute constant offset
	*objoffset = 0;
	for(i = 0; i < num_firststgcols; i++) *objoffset = *objoffset + *(kappa + i) * (*(u0 + i)) * (*(u0 + i));
	*objoffset = .5*mu * *objoffset;
	
	// compute linear terms
	if(simplifybinaries == 1){
		for(i = 0; i < num_firststgcols; i++){
			if( *(coltype + i) == 'B' ){
				*(obj + i) = *(msumx + i) - mu * (*(u0 + i) - .5);	// (eta = 1, kappa = 1, see above)
			} else {
				*(obj + i) = (*(eta + i)) * (*(msumx + i)) - mu * *(kappa + i) * *(u0 + i);
			}
		}
	} else {
		for(i = 0; i < num_firststgcols; i++)
			*(obj + i) = (*(eta + i)) * (*(msumx + i)) - mu * *(kappa + i) * *(u0 + i);
		// I don't separe binaries from non binary variables here to take advantage of optimized caching
	}
	
	// compute quadratic coefficients
	if(simplifybinaries == 1){
		j = 0;
		for(i = 0; i < ncolsxprs && j < nqtr; i++){
			if( *(coltype + i) != 'B' ){
				*(mqc1 + j) = i;
				*(mqc2 + j) = i;
				*(dqe + j) = mu * *(kappa + i);
				j++;
			}
		}
	} else {
		for(i = 0; i < nqtr; i++){
			*(mqc1 + i) = i;
			*(mqc2 + i) = i;
			*(dqe + i) = mu * *(kappa + i);
		}
	}
	
	// load quadratic objective into the optimizer
	int *clind = (int*) malloc(sizeof(int)*ncolsxprs);
	for(i = 0; i < ncolsxprs; i++) *(clind + i) = i;
	XPRSchgobj(*firststageprob, ncolsxprs, clind, obj);
	XPRSchgmqobj(*firststageprob, nqtr, mqc1, mqc2, dqe);
	
	// free heap
	free(coltype);
	free(eta);
	free(kappa);
	free(obj);
	free(mqc1);
	free(mqc2);
	free(dqe);
	free(clind);
	
	// return success indicator
	return(0);
}

/* Function to load objective for non-scenario MILP subproblem: f_0(-sum x) = min_{u in U} {((-sum x)./(u_UB - u_LB))^T u} */
int load_f0_lin_objective(const int num_firststgcols, const double *msumx,
	const double *firststagelb, const double *firststageub,
	XPRSprob *firststageprob)
{
	// iterators
	int i;
	
	// check that problem dimensions agree with firts stage variables
	int ncolsxprs;
	XPRSgetintattrib(*firststageprob, XPRS_ORIGINALCOLS, &ncolsxprs);
	if( ncolsxprs != num_firststgcols ){
		printf("\nProvided XPRSprob (firststageprob) has a number of columns (%d) incompatible"
			"\nwith the provided number of first stage columns (%d).\n\n",
			ncolsxprs, num_firststgcols);
		return(-1);
	}
	
	// load new linear terms
	int *clind = (int*) malloc(sizeof(int)*ncolsxprs);
	for(i = 0; i < ncolsxprs; i++) *(clind + i) = i;
	double *msumxscaled = (double*) malloc(sizeof(double)*ncolsxprs);
	for(i = 0; i < ncolsxprs; i++){
		if( (*(firststageub + i)) - (*(firststagelb + i)) > DBL_EPSILON ){
			*(msumxscaled + i) = 1.0/((*(firststageub + i)) - (*(firststagelb + i))) * (*(msumx + i));
		} else {
			*(msumxscaled + i) = *(msumx + i);
		}
	}
	XPRSchgobj(*firststageprob, ncolsxprs, clind, msumxscaled);
	free(clind);
	free(msumxscaled);
	
	// return success indicator
	return(0);
}

/* Function to erase current objective of first stage problem */
int clean_f0_objective(XPRSprob *firststageprob)
{
	// iterators
	int i, j;
	
	// get number of variables
	int ncolsxprs;
	XPRSgetintattrib(*firststageprob, XPRS_ORIGINALCOLS, &ncolsxprs);
	
	// delete current quadratic terms if they exist
	int nqelems;
	XPRSgetmqobj(*firststageprob, NULL, NULL, NULL, 0, &nqelems, 0, ncolsxprs-1);
	if(nqelems > 0){
		int nqelemsret, *mstart, *mrwind, *mclind;
		double *dmatval;
		mstart = (int*) malloc(sizeof(int)*(ncolsxprs+1));
		mclind = (int*) malloc(sizeof(int)*nqelems);
		XPRSgetmqobj(*firststageprob, mstart, mclind, NULL, nqelems, &nqelemsret, 0, ncolsxprs-1);
		if( nqelems != nqelemsret || *(mstart + ncolsxprs) != nqelemsret ){
			printf("\nXpress is returning an erratic number of elements in the quadratic matrix.\n");
			return(-1);
		}
		mrwind = (int*) malloc(sizeof(int)*nqelems);
		dmatval = (double*) malloc(sizeof(double)*nqelems);
		for(i = 0; i < ncolsxprs; i++){
			for(j = *(mstart + i); j < *(mstart + i + 1); j++){
				*(mrwind + j) = i;
				*(dmatval + j) = 0.0;
			}
		}
		XPRSchgmqobj(*firststageprob, nqelems, mrwind, mclind, dmatval);
		free(mstart);
		free(mclind);
		free(mrwind);
		free(dmatval);
	}
	
	// erase current linear terms
	int *clind = (int*) malloc(sizeof(int)*ncolsxprs);
	double *obj = (double*) malloc(sizeof(double)*ncolsxprs);
	for(i = 0; i < ncolsxprs; i++) *(clind + i) = i;
	for(i = 0; i < ncolsxprs; i++) *(obj + i) = 0;
	XPRSchgobj(*firststageprob, ncolsxprs, clind, obj);
	free(clind);
	free(obj);
	
	// return success indicator
	return(0);
}
