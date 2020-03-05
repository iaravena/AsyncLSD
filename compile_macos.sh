#!/bin/bash

# get into source directory
cd source

# compile and link
mpicc -c Asynchronous.c Coordinator.c Coordinator_DualJobs.c Coordinator_PrimalJobs.c Coordinator_MetaProc.c Coordinator_Misc.c Worker.c HighLevelFunc.c LowLevelFunc.c ReadOptions.c SMPS_prepare.c SMPS_Xpress_Read.c SMPS_Xpress_ScenSubprob.c SMPS_Xpress_StageSubprobs.c -O2 -g -D_REENTRANT -I${XPRESSDIR}/include
mpicc Asynchronous.o Coordinator.o Coordinator_DualJobs.o Coordinator_PrimalJobs.o Coordinator_MetaProc.o Coordinator_Misc.o Worker.o HighLevelFunc.o LowLevelFunc.o ReadOptions.o SMPS_prepare.o SMPS_Xpress_Read.o SMPS_Xpress_ScenSubprob.o SMPS_Xpress_StageSubprobs.o -o AsyncLSD -L${XPRESSDIR}/lib -lxprs -lm -g

# delete intermediate files
rm Asynchronous.o Coordinator.o Coordinator_DualJobs.o Coordinator_PrimalJobs.o Coordinator_MetaProc.o Coordinator_Misc.o Worker.o HighLevelFunc.o LowLevelFunc.o ReadOptions.o SMPS_prepare.o SMPS_Xpress_Read.o SMPS_Xpress_ScenSubprob.o SMPS_Xpress_StageSubprobs.o

# make executable
chmod +x AsyncLSD

# go back to the original directory
mv AsyncLSD ../
cd ..
