REM get into source directory
cd source

REM compile and link
gcc -c Asynchronous.c Coordinator.c Coordinator_DualJobs.c Coordinator_PrimalJobs.c Coordinator_MetaProc.c Coordinator_Misc.c Worker.c HighLevelFunc.c LowLevelFunc.c ReadOptions.c SMPS_prepare.c SMPS_Xpress_Read.c SMPS_Xpress_ScenSubprob.c SMPS_Xpress_StageSubprobs.c -O2 -I%MSMPI_INC_GCC%/ -D_REENTRANT -I%XPRESSDIR%/include
gcc Asynchronous.o Coordinator.o Coordinator_DualJobs.o Coordinator_PrimalJobs.o Coordinator_MetaProc.o Coordinator_Misc.o Worker.o HighLevelFunc.o LowLevelFunc.o ReadOptions.o SMPS_prepare.o SMPS_Xpress_Read.o SMPS_Xpress_ScenSubprob.o SMPS_Xpress_StageSubprobs.o -o AsyncLSD.exe -L%MSMPI_LIB64_GCC%/ -lmsmpi -L%XPRESSDIR%/lib -lxprs

REM delete intermediate files
del Asynchronous.o Coordinator.o Coordinator_DualJobs.o Coordinator_PrimalJobs.o Coordinator_MetaProc.o Coordinator_Misc.o Worker.o HighLevelFunc.o LowLevelFunc.o ReadOptions.o SMPS_prepare.o SMPS_Xpress_Read.o SMPS_Xpress_ScenSubprob.o SMPS_Xpress_StageSubprobs.o

REM go back to the original directory
move AsyncLSD.exe ../
cd ..
