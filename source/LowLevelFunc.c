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
 * Low level functions: perform basic operations such as sorting, comparison, etc...
 */

/* Low level header */
#include "LowLevelFunc.h"

/* OSX headers for accurate time measuring */
#if (defined __APPLE__ || defined __MACH__)
	#include <mach/clock.h>
	#include <mach/mach.h>
#endif

/* Function to get index of a string within a buffer */
int get_string_index(const struct string_buffer *sbuffer, const char *test_string)
{
	int i;
	char test_string_filled[(*sbuffer).atom_len+1];
	
	// trivial not found cases
	if((*sbuffer).num_elems == 0 || strlen(test_string) > (*sbuffer).atom_len)
		return(-1);
		
	// add blank spaces if necessary
	add_spaces(test_string, (*sbuffer).atom_len+1, (*sbuffer).atom_len, &test_string_filled[0]);
	
	// loop over the buffer, looking for the test string
	i = 0;
	while(i < (*sbuffer).num_elems){
		if(strcmp(test_string_filled, (*sbuffer).buffer + ((*sbuffer).atom_len+1)*i) == 0){
			break;
		}
		i++;
	}
	
	if( i < (*sbuffer).num_elems ){
		return(i);
	} else {
		return(-1);
	}
}

/* Function to insert a string on a predefined position in a buffer */
int insert_in_buffer(struct string_buffer *sbuffer, const int position, const char *new_atom)
{
	int len = strlen(new_atom);
	char new_atom_filled[(*sbuffer).atom_len+1];
	
	// inconsistent insertion
	if(position >= (*sbuffer).num_elems || len > (*sbuffer).atom_len)
		return(-1);
	
	// add blank spaces if necessary
	add_spaces(new_atom, (*sbuffer).atom_len+1, (*sbuffer).atom_len, &new_atom_filled[0]);
	
	// inset new atom in the desired position
	strcpy((*sbuffer).buffer+((*sbuffer).atom_len+1)*position, new_atom_filled);
	
	// return success indicator
	return(0);
}

/* Function to get index of a string within a buffer, inserting the string at the tail if it is not found */
int get_index_within_buffer(struct string_buffer *sbuffer, const char *test_atom)
{
	int i;
	int len = strlen(test_atom);
	char test_atom_filled[(*sbuffer).atom_len+1];
	
	// inconsistent insertion
	if(len > (*sbuffer).atom_len)
		return(-1);
	
	// add blank spaces if necessary
	add_spaces(test_atom, (*sbuffer).atom_len+1, (*sbuffer).atom_len, &test_atom_filled[0]);
	
	i = 0;
	while(i < (*sbuffer).num_elems){
		if(strcmp(test_atom_filled, (*sbuffer).buffer + ((*sbuffer).atom_len+1)*i) == 0){
			break;
		}
		i++;
	}
	
	// if atom is not found, insert it at the tail
	if(i == (*sbuffer).num_elems){
		strcpy((*sbuffer).buffer+((*sbuffer).atom_len+1)*i, test_atom_filled);
		(*sbuffer).num_elems = i+1;
	}
	
	// return success indicator
	return(i);
}

/* Function to add trailing spaces to a string
 * source: http://stackoverflow.com/questions/21855807/c-how-to-append-concatenate-x-spaces-to-a-string */
void add_spaces(const char *source, const int size, int num_of_spaces, char *dest)
{
	int len = strlen(source);
	strcpy(dest, source);
	if( len + num_of_spaces >= size ) {
		num_of_spaces = size - len - 1;
	}
	memset(dest+len, ' ', num_of_spaces);
	dest[len + num_of_spaces] = '\0';
}

/* Function to add trailing spaces to a string in-place */
char* add_spaces_in_place(char *strptr, const int size, int num_of_spaces)
{
	int len = strlen(strptr);
	if( len + num_of_spaces >= size ) {
		num_of_spaces = size - len - 1;
	}
	memset(strptr+len, ' ', num_of_spaces);
	strptr[len + num_of_spaces] = '\0';
	return(strptr);
}

/* Function to find next non-empty character position in a string */
int next_nonempty_pos(char *strptr, const int size, int position)
{
	do{
		position++;
	} while (*(strptr + position)==' ' && position < size - 1);
	return(position);
}

/* Function to see whether the elements of an integer vector are contained into a larger vector
 * lookfor and lookin are assumed to be sorted in increasing order, and have only unique elements */
int integer_set_inclusion(const size_t lenlookfor, const int *lookfor,
	const size_t lenlookin, const int *lookin)
{
	// check if the lenghts can let us out of here quickly
	if(lenlookfor > lenlookin){
		return(0);
	}
	
	// check element by element of lookfor
	int i, j = 0;
	for(i = 0; i < lenlookfor; i++){
		// increase j until find a hit
		while( *(lookfor + i) > *(lookin + j) && j < lenlookin ) j++;
		// check whether match is exact, if not we can return
		if( *(lookfor + i) != *(lookin + j) ){
			return(0);
		}
	}
	
	// if we get here -> all elements in lookfor have been matched exactly by an element in lookin
	return(1);
}

/* Function to compare integers */
int compare_integers(const void *a, const void *b){
	return ( *((int*) a) - *((int*) b) );
}

/* Function to compare integers in a vector */
int compare_integers_r(const void *i1_arg, const void *i2_arg, void *vector_arg){
	int i1 = *((int*) i1_arg);
	int i2 = *((int*) i2_arg);
	int *vector = (int*) vector_arg; 
	return ( (*(vector + i1)) - (*(vector + i2)) );
}

/* Compare two strings stored in a string buffer */
int compare_atoms_buffer(const void *i1_arg, const void *i2_arg, void *sbuffer_arg){
	int i1 = *((int*) i1_arg);
	int i2 = *((int*) i2_arg);
	struct string_buffer *sbuffer = (struct string_buffer *) sbuffer_arg;
	return(strcmp((*sbuffer).buffer + ((*sbuffer).atom_len+1) * i1,
		(*sbuffer).buffer + ((*sbuffer).atom_len+1) * i2));
}

/* Lexicographically compare rows of an integer data frame */
int compare_rows_idf(const void *i1_arg, const void *i2_arg, void *df_arg){
	int i1 = *((int*) i1_arg);
	int i2 = *((int*) i2_arg);
	struct data_frame *df = (struct data_frame *) df_arg;
	int j = 0, diff = 0;
	do {
		diff = *(*((int **) (*df).ptr + j) + i1) - *(*((int **) (*df).ptr + j) + i2);
		j++;
	} while (j < (*df).numcols && diff == 0);
	return(diff);
}

/* Multi-platform quick sort: dispatchs either qsort_s or qsort_r
 * source: http://stackoverflow.com/questions/4300896/how-portable-is-the-re-entrant-qsort-r-function-compared-to-qsort */
struct sort_r_data
{
	void *arg;
	int (*compar)(const void *a1, const void *a2, void *aarg);
} sort_r_data;
int sort_r_arg_swap(void *s, const void *aa, const void *bb)
{
	struct sort_r_data *ss = (struct sort_r_data*)s;
	return (ss->compar)(aa, bb, ss->arg);
}
void sort_r(void *base, size_t nel, size_t width,
	int (*compar)(const void *a1, const void *a2, void *aarg), void *arg)
{
	#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux__)
		qsort_r(base, nel, width, compar, arg);
	#elif (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || \
		defined __FREEBSD__ || defined __BSD__ || \
		defined OpenBSD3_1 || defined OpenBSD3_9)
		struct sort_r_data tmp;
		tmp.arg = arg;
		tmp.compar = compar;
		qsort_r(base, nel, width, &tmp, &sort_r_arg_swap);
	#elif (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		struct sort_r_data tmp = {arg, compar};
		qsort_s(base, nel, width, &sort_r_arg_swap, &tmp);
	#else
		#error Cannot detect operating system
	#endif
}

/* Skipping lines when reading a file */
int flineskip(FILE *filecon, const int nlines, const int buffersize){
	int i = 0;
	char linetoignore[buffersize];
	// skip nlines
	while(i < nlines
		&& fgets(linetoignore,sizeof linetoignore, filecon) !=  NULL)
		i++;
	// return difference between number of lines asked to skip and skipped lines
	if(i < nlines){
		return(nlines - (i+1));
	} else {
		return(0);
	}
}

/* Get indices of elements of a buffer within another, presorted, buffer (sorted with compare_atoms_buffer)*/
void match_to_presorted_buffer(const struct string_buffer *lookin_buffer,
	const int *linbuff_sorted_ix, const struct string_buffer *lookfor_buffer,
	int *result_indices){
	
	int i, j;
	int *lfrbuff_sorted_ix;
	char lookforstr[(*lookin_buffer).atom_len+1];
	
	// generate array of sorted look for elements
	lfrbuff_sorted_ix = (int *) calloc((*lookfor_buffer).num_elems, sizeof(int));
	for(i = 0; i < (*lookfor_buffer).num_elems; i++) *(lfrbuff_sorted_ix + i) = i;
	sort_r(lfrbuff_sorted_ix, (*lookfor_buffer).num_elems, sizeof(int),
		&compare_atoms_buffer, (void *) lookfor_buffer);
	
	// find elements of lookfor in lookin
	j = 0;
	for(i = 0; i < (*lookfor_buffer).num_elems; i++){
		// format the string we look for
		//printf("%d\t%s\n", i, (*lookfor_buffer).buffer + ((*lookfor_buffer).atom_len+1) * *(lfrbuff_sorted_ix+i));
		add_spaces((*lookfor_buffer).buffer + ((*lookfor_buffer).atom_len+1) * *(lfrbuff_sorted_ix+i),
			(*lookin_buffer).atom_len+1, (*lookin_buffer).atom_len, &lookforstr[0]);
		// search lookin incrementally, until the end
		while(j < (*lookin_buffer).num_elems && strcmp(lookforstr,
			(*lookin_buffer).buffer + ((*lookin_buffer).atom_len+1) * *(linbuff_sorted_ix + j)) > 0){
			j++;
		}
		// assign indices
		if(j < (*lookin_buffer).num_elems && strcmp(lookforstr,
			(*lookin_buffer).buffer + ((*lookin_buffer).atom_len+1) * *(linbuff_sorted_ix + j)) == 0){
			*(result_indices + *(lfrbuff_sorted_ix+i)) = *(linbuff_sorted_ix + j);
		} else {
			*(result_indices + *(lfrbuff_sorted_ix+i)) = -1;
		}
	}
	
	// free memory
	free(lfrbuff_sorted_ix);
	
}

/* Swap two elements in an integer vector */
void swap_elements_int_vector(int *vector, const int i, const int j)
{
	int valuei = *(vector + i);
	*(vector + i) = *(vector + j);
	*(vector + j) = valuei;
}

/* Reorder vector using given indices (indices will be render useless) */
void reorder_int_vector(const size_t numelems, int *unsrtvec, int *indices)
{
	// iterators
	int i, j;
	// sort rows
	for (i = 0; i < numelems; i++){
		// swap rows at i and *(indices + i)
		swap_elements_int_vector(unsrtvec, i, *(indices + i));
		// now, tell the sorted indices the location of the old item after swap
		for (j = i; j < numelems; j++){
			if( *(indices + j) == i ){
				*(indices + j) = *(indices + i);
				break; // we only need the first one, so then we're done
			}
		}
	}
}

/* Swap rows for a column of a data frame */
void swap_rows_incol_data_frame(struct data_frame *df, const int col, const int rowi, const int rowj)
{
	char *charcol, charelem;
	int *intcol, intelem;
	double *dblcol, dblelem;
	struct string_buffer *strbuffcol;
	char *strbuffelem;
	
	switch(*((*df).coltypes + col)){
		case DFR_COLTYPE_CHAR :
			charcol = (char *) *((*df).ptr + col);
			charelem = *(charcol + rowi);
			*(charcol + rowi) = *(charcol + rowj);
			*(charcol + rowj) = charelem;
			break;
		case DFR_COLTYPE_INT :
			intcol = (int *) *((*df).ptr + col);
			intelem = *(intcol + rowi);
			*(intcol + rowi) = *(intcol + rowj);
			*(intcol + rowj) = intelem;
			break;
		case DFR_COLTYPE_DOUBLE :
			dblcol = (double *) *((*df).ptr + col);
			dblelem = *(dblcol + rowi);
			*(dblcol + rowi) = *(dblcol + rowj);
			*(dblcol + rowj) = dblelem;
			break;
		case DFR_COLTYPE_STRBUFF :
			strbuffcol = (struct string_buffer*) *((*df).ptr + col);
			strbuffelem = (char *) malloc(sizeof(char)*((*strbuffcol).atom_len+1));
			strcpy(strbuffelem, (*strbuffcol).buffer + ((*strbuffcol).atom_len+1)*rowi);
			strcpy((*strbuffcol).buffer + ((*strbuffcol).atom_len+1)*rowi,
				(*strbuffcol).buffer + ((*strbuffcol).atom_len+1)*rowj);
			strcpy((*strbuffcol).buffer + ((*strbuffcol).atom_len+1)*rowj, strbuffelem);
			free(strbuffelem);
			break;
	}
}

/* Swap two rows of a data frame */
void swap_rows_data_frame(struct data_frame *df, const int rowi, const int rowj)
{
	int j;
	if(rowi != rowj){
		for(j = 0; j < (*df).numcols; j++) swap_rows_incol_data_frame(df, j, rowi, rowj);
	}
}

/* Reorder data frame using given indices (indices will be render useless) */
void reorder_data_frame(struct data_frame *unsrtdf, int *indices)
{
	// iterators
	int i, j, newIndex;
	// sort rows
	for (i = 0; i < (*unsrtdf).numrows; i++){
		// store the current index
		newIndex = *(indices + i);
		// swap rows at newIndex and i
		swap_rows_data_frame(unsrtdf, newIndex, i);
		// now, tell the sorted indices the location of the old item after swap
		for (j = i; j < (*unsrtdf).numrows; j++){
			if (indices[j] == i){
				indices[j] = newIndex;
				break; // we only need the first one, so then we're done
			}
		}
	}
}

/* Function to get row or column indices associated with a set of time stages */
void get_items_by_property_value(const size_t numitems, const int *items,
	const int *itempropertyarr, const size_t numadmvalues, int *admvalues,
	const int maxreturnitems, int *numfounditems, int *founditems)
{
	// iterators
	int i, j, k;
	
	// sort admissible values (in place)
	qsort((void*) admvalues, numadmvalues, sizeof(int), &compare_integers);
	
	// sort property array
	struct data_frame propdf;
	void *propcol = (void*) itempropertyarr;
	propdf.numcols = 1;
	propdf.numrows = numitems;
	propdf.ptr = &propcol;
	int *itempropsort = malloc(sizeof(int)*numitems);
	for(i = 0; i < numitems; i++) *(itempropsort + i) = i;
	sort_r((void*) itempropsort, numitems, sizeof(int), &compare_rows_idf, (void*) &propdf);
	
	// find items (if maxreturnitems is 0, just count the number of items)
	int delta;
	i = 0;
	j = 0;
	k = 0;
	for(i = 0; i < numitems; i++){
		delta = *(itempropertyarr + *(itempropsort + i)) - *(admvalues + j);
		// increase counter j until admissible value reachs current property values
		while(delta > 0 && j < numadmvalues - 1){
			j++;
			delta = *(itempropertyarr + *(itempropsort + i)) - *(admvalues + j);
		}
		// check wether we have found a new item and record it
		if(delta == 0){
			if(k < maxreturnitems){
				*(founditems + k) = *(items + *(itempropsort + i));
			}
			k++;
		}
		// check whether we should continue looking
		if(delta > 0 && j == numadmvalues - 1) break;
	}
	*numfounditems = k;
	
	// sort found items
	if(maxreturnitems > 0){
		qsort((void*) founditems, (int) fmin(maxreturnitems, *numfounditems), sizeof(int),
			&compare_integers);
	}
	
	// free heap
	free(itempropsort);
}

/* Work around lack of clock_gettime in OS X
 * source: https://gist.github.com/jbenet/1087739 */
void current_utc_time(struct timespec *ts)
{
	#if (defined __APPLE__ || defined __MACH__) // OS X does not have clock_gettime, use clock_get_time
	  clock_serv_t cclock;
	  mach_timespec_t mts;
	  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	  clock_get_time(cclock, &mts);
	  mach_port_deallocate(mach_task_self(), cclock);
	  ts->tv_sec = mts.tv_sec;
	  ts->tv_nsec = mts.tv_nsec;
	#else
	  clock_gettime(CLOCK_REALTIME, ts);
	#endif
}

/* Random (approximatedly) Bernoulli number generator */
int rbernoulli(const double succesprob)
{
	return(((double) rand())/RAND_MAX <= succesprob);
}

/* Random sample with probabilities */
void nonuniformsample(const int numelems, const double *probabilities,
	const int samplesize, int *sampleindices, int *frequency, int *tablesize)
{
	// compute the cummulative probabilites
	double *cumprobabilites = (double*) malloc(sizeof(double)*numelems);
	{
		int i;
		*cumprobabilites = *probabilities;
		for(i = 1; i < numelems; i++){
			*(cumprobabilites + i) = *(cumprobabilites + i - 1) + *(probabilities + i);
		}
	}
	
	// scale cumulative probabilities to 1
	if( *(cumprobabilites + numelems - 1) != 1 ){
		int i;
		double probscale = 1.0 / *(cumprobabilites + numelems - 1);
		for(i = 0; i < numelems; i++){
			*(cumprobabilites + i) *= probscale;
		}
	}
	
	// output format: long vector or frequency table
	int *internalindices;
	if( frequency == NULL || tablesize == NULL ){
		// long vector
		internalindices = sampleindices;
	}else{
		// frequency table
		internalindices = (int*) malloc(sizeof(int)*samplesize);
		*tablesize = 0;
	}
	
	// sample elements and sort them
	{
		int i;
		for(i = 0; i < samplesize; i++){
			*(internalindices + i) = rand();
		}
		qsort(internalindices, samplesize, sizeof(int), &compare_integers);
	}
	
	// find break points and assing true indices
	{
		int i, j, jstart = 0, jend = 0;
		for(i = 0; i < numelems; i++){
			while( ((double) *(internalindices + jend)) <= RAND_MAX * *(cumprobabilites + i)
				&& jend < samplesize ){
				jend++;
			}
			if( jstart < jend ){
				if( frequency == NULL || tablesize == NULL ){
					// long vector
					for(j = jstart; j < jend; j++){
						*(sampleindices + j) = i;
					}
				}else{
					// frequency table
					*(sampleindices + *tablesize) = i;
					*(frequency + *tablesize) = jend - jstart;
					*tablesize += 1;
				}
				jstart = jend;
			}
		}
	}
	
	// release allocated memory
	free(cumprobabilites);
	if( !(frequency == NULL || tablesize == NULL) ){
		free(internalindices);
	}
}

/* L1 norm */
double L1norm(size_t dimension, const double *vec1, const double *vec2)
{
	size_t i;
	double normvalue = 0;
	for(i = 0; i < dimension; i++){
		normvalue += fabs( *(vec1 + i) - *(vec2 + i) );
	}
	return(normvalue);
}

/* Cycling next element: get next element from the ring {0,...,ringsize-1} */
int cycling_next_element(const int ringsize, const int current_element)
{
	// Check that input makes sense
	if(current_element < 0 || current_element >= ringsize){
		printf("\nCurrent element is not in the ring.\n");
		return(-1);
	}
	
	// get next element
	if( current_element + 1 < ringsize ){
		return(current_element + 1);
	}else{
		return(0);
	}
}

/* Cycling difference: compute difference between numbers in the ring {0,...,ringsize-1}
int cycling_difference(const int ringsize, const int substractfrom, const int substractwhat)
{
	// Check that input makes sense
	if(substractfrom < 0 || substractfrom >= ringsize
		|| substractwhat < 0 || substractwhat >= ringsize){
		printf("\nElements to substract are not in the ring.\n");
		return(-1);
	}
	
	// Compute difference
	if(substractwhat <= substractfrom){		// return ordinary difference
		return(substractfrom - substractwhat);
	}else{									// return difference over the ring
		return(ringsize - substractwhat + substractfrom);
	}
	
}
*/
