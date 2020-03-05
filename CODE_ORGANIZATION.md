# Source code organization

The source code is organized in a layered structure, from low-level operations like sorting a character buffer to high-level operations like solving a subproblem for a certain scenario. In the following we briefly describe the source code from the higher to the lower level.

**NOTE**: The source code is commented step-by-step, for further clarifications take a look at the source code itself.

## Top level

- [Asynchronous.c](source/Asynchronous.c) : Contains the 'main' function, which launches either the Coordinator (*Master*) or Worker (*Slave*) procedure depending on the MPI rank of the process.
- [Coordinator.c](source/Coordinator.c) : Implements the *Master*, presented schematically in Figure 4 of [the paper](README.md#citation). Subroutines of the *Master*, almost always corresponding to boxes in Figure 4, have been spread across four files for readability;
	* [Coordinator_MetaProc.c](source/Coordinator_MetaProc.c) : Contains routines to allocate memory and executing the Initialization phase (section 5.3 of [the paper](README.md#citation)).
	* [Coordinator_DualJobs.c](source/Coordinator_DualJobs.c) : Contains routines to launch and post-process dual tasks, 'dual f0' and 'dual scenario' in Figure 4.
	* [Coordinator_PrimalJobs.c](source/Coordinator_PrimalJobs.c) : Contains routines to launch and post-process primal tasks, 'second stage scenario' and 'primal projection' in Figure 4.
	* [Coordinator_Misc.c](source/Coordinator_Misc.c) : Contains functions used by the *Master* to determine the next task, compute the current stepsize, among others.
- [Worker.c](source/Worker.c) : Implements the *Slave*, presented schematically in Figure 3 of [the paper](README.md#citation).

## High level
	
- [HighLevelFunc.c](source/HighLevelFunc.c) : High level functions for reading SMPS instances, transforming the CORE problem (held in memory by each *Slave*) into scenario subproblems and vice-versa, and defining MPI message types. These functions, and functions at the [Top level](#top-level), include calls to `malloc` that are not released within the same function, which is not the case for functions at lower levels.

## Mid level
	
- [ReadOptions.c](source/ReadOptions.c) : Implements a function for reading the configuration parameters file.
- [SMPS_Xpress_Read.c](source/SMPS_Xpress_Read.c) : Implements functions for reading SMPS TIME and STOCH files.
- [SMPS_Xpress_ScenSubProb.c](source/SMPS_Xpress_ScenSubProb.c) : Implements functions for modifying CORE or scenario problems held in memory by Xpress, including changing objective and constraints, fixing variables and performing the period relaxation (used during Initialization).
- [SMPS_Xpress_StageSubProbs.c](source/SMPS_Xpress_StageSubProbs.c) : Implements functions for creating first stage problems, necessary for evaluating f0 and f0^mu, and period subproblems, from CORE or scenario subproblems held in memory by Xpress.
	
## Low level
	
- [LowLevelFunc.c](source/LowLevelFunc.c) : Implements functions to manage character buffers, manage 'data frames' (tables stored in a single struct), keeping track of time in different operating systems, sorting in place different operating systems, among others.
