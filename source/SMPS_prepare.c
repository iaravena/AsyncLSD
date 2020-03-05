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
 * Function for creating folders, decompress instance and prepare files for execution
 */

/* SMPS Xpress algorithm header */
#include "SMPS_Xpress.h"

/* Function to create execution folder and prepare files */
int prepare_files(const char *workdir, const char *instdir, const char *SMPSrootname,
	const char *optionsfile, const struct options *opt)
{
	// Define variables
	char syscmd[1024];
	char CORfname[200], MPSfname[200];
	
	// Create working directory
	#if (defined _WIN32 || defined _WIN64 || defined __Wload_INDOWS__)
	// Windows code
		mkdir(workdir);
	#else
	// GNU/Linux code
		mkdir(workdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	#endif
	
	// Decompress SMPS files in the working directory
	// Windows code (7zip commands, just to debug)
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
		// decompress tar.gz, write tar to workdir
		strcpy(syscmd, WIN_7ZIP_EXECUTABLE " e -o");
		strcat(strcat(syscmd, workdir), " ");
		strcat(strcat(strcat(strcat(syscmd, instdir), SYSTEM_SLASH), SMPSrootname), ".tar.gz > nul");
		system(syscmd);
		// extract tar contents
		syscmd[0] = '\0';
		strcpy(syscmd, WIN_7ZIP_EXECUTABLE " e -o");
		strcat(strcat(syscmd, workdir), " ");
		strcat(strcat(strcat(strcat(syscmd, workdir), SYSTEM_SLASH), SMPSrootname), ".tar > nul");
		system(syscmd);
		// delete tar ball
		syscmd[0] = '\0';
		strcpy(syscmd, "del ");
		printf("\n%s\n", strcat(strcat(strcat(strcat(syscmd, workdir), SYSTEM_SLASH), SMPSrootname), ".tar"));
		system(syscmd);
	#else
	// GNU/Linux code
		// extract tar.gz contents to workdir
		strcpy(syscmd, "tar -xf ");
		strcat(strcat(strcat(strcat(syscmd, instdir), SYSTEM_SLASH), SMPSrootname), ".tar.gz");
		strcat(strcat(syscmd, " -C "), workdir);
		system(syscmd);
	#endif
	
	// Rename CORE file as MPS
	strcpy(CORfname, workdir);
	strcat(strcat(strcat(CORfname, SYSTEM_SLASH), SMPSrootname), ".cor");
	strcpy(MPSfname, workdir);
	strcat(strcat(strcat(MPSfname, SYSTEM_SLASH), SMPSrootname), ".mps");
	rename(CORfname, MPSfname);
	
	// Copy options file to the working directory
	#if (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)	// Windows code
		strcpy(syscmd, "copy ");
	#else															// GNU/Linux code
		strcpy(syscmd, "cp ");
	#endif
	strcat(strcat(syscmd, optionsfile), " ");
	strcat(strcat(syscmd, workdir), SYSTEM_SLASH "execution_options.txt");
	system(syscmd);
	
	/* Return success indicator */
	return(0);
}