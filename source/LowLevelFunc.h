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

#ifndef LOW_LEVEL_INCLUDED
	
	// define LOW_LEVEL_INCLUDED to avoid loading the header again
	#define LOW_LEVEL_INCLUDED
	
	/* INCLUDE STATEMENTS */
	
	// Standard libraries
	#include <stdio.h>			// provides functions for reading and printing
	#include <stdlib.h>			// provides several functions for managing memory, files and interact with system
	#include <string.h>			// managing strings
	#include <math.h>			// math functions
	#include <time.h>			// time functions
	
	/* CONSTANTS */
	
	// Dataframe column types
	#define DFR_COLTYPE_CHAR			101
	#define DFR_COLTYPE_INT				102
	#define DFR_COLTYPE_DOUBLE			103
	#define DFR_COLTYPE_STRBUFF			104		// string_buffer structure
	
	// Time measuring
	#define NANOSEC2SEC					1E-9
	
	/* STRUCTS */
	typedef struct string_buffer
	{
		int atom_len, num_elems;
		char *buffer;
	} string_buffer;
	typedef struct data_frame
	{
		int numcols, numrows;
		void **ptr;
		int *coltypes;
	} data_frame;
	
	/* FUNCTIONS AND PROCEDURES */
	extern int get_string_index(const struct string_buffer*, const char*);
	extern void add_spaces(const char*, const int, int, char*);
	extern char* add_spaces_in_place(char*, const int, int);
	extern int insert_in_buffer(struct string_buffer*, const int, const char*);
	extern int get_index_within_buffer(struct string_buffer*, const char*);
	extern int next_nonempty_pos(char*, const int, int);
	extern int integer_set_inclusion(const size_t, const int*, const size_t, const int*);
	extern int compare_integers(const void*, const void*);
	extern int compare_integers_r(const void*, const void*, void*);
	extern int compare_atoms_buffer(const void*, const void*, void*);
	extern int compare_rows_idf(const void*, const void*, void*);
	extern void sort_r(void*, size_t, size_t, int (*)(const void*, const void*, void*), void*);
	extern int flineskip(FILE*, const int, const int);
	extern void match_to_presorted_buffer(const struct string_buffer*, const int*,
		const struct string_buffer*, int*);
	extern void reorder_int_vector(const size_t, int*, int*);
	//extern void swap_rows_data_frame(struct data_frame*, const int, const int);
	extern void reorder_data_frame(struct data_frame*, int*);
	extern void get_items_by_property_value(const size_t, const int*,
		const int*, const size_t, int*, const int, int*, int*);
	extern void current_utc_time(struct timespec*);
	extern int rbernoulli(const double);
	extern void nonuniformsample(const int, const double*, const int, int*, int*, int*);
	extern double L1norm(size_t, const double*, const double*);
	extern int cycling_next_element(const int, const int);
	//extern int cycling_difference(const int, const int, const int);
	
#endif