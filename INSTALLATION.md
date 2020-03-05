# Installation instructions

## Prerrequisites

### Linux/Mac OS

You will need the following programs installed on your system before you can compile the code or execute the asynchronous algorithm:

- **C compiler**: We recommend a recent version of the GNU Compiler Collection (GCC 4.9 or superior). It comes included in most Linux distributions. Mac OS users can obtain GCC as part of the 'Command Line Tools for Xcode'. Other recent C compilers might also be used.
- **MPI implementation**: We recommend a recent version of Open MPI (1.10.5 or superior). The software can be obtained for free at <https://www.open-mpi.org/>. Other implementations of MPI can also be used as long as they support the MPI 3.0 Standard.
  *IMPORTANT*: Once you are done with the installation, make sure that you can compile MPI code using the `mpicc` command.
- **FICO Xpress 7.8** or superior: Xpress is a commercial solver for mathematical programs. You will need an Xpress license before you can install or download the solver. If you are an academic user, you can request an academic license (free of charge) at <http://subscribe.fico.com/Academic-Partner-Program>. Once you have downloaded and installed Xpress, append the following lines to your `.bashrc` or `.bash_profile` file:
```
export XPRESSDIR=/path/to/xpressmp
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XPRESSDIR}/lib
export PATH=${PATH}:${XPRESSDIR}/bin
```

### Windows

- **MinGW 6.2.0** or superior: MinGW is recommended here as is one of the simplest way to obtain GCC Windows. The software can be obtained for free from <http://www.mingw.org/>. Version 6.2.0 or superior is required by Microsoft MPI. Once you have downloaded and placed MinGW in its definitive folder, locate the bin directory within your MinGW installation folder. Append this directory to the `%Path%` environment variable (Control Panel > System > Advanced system settings > Advanced > Environment Variables).
- **Microsoft MPI 7.1** or superior: Microsoft MPI (MSMPI) is a Microsoft implementation of the Message Passing Interface standard for developing and running parallel applications on the Windows platform. It can be obtained for free from <https://msdn.microsoft.com/en-us/library/bb524831>. Follow the instructions for download and installation.
	In order to use MSMPI with GCC we will need to set additional environment variables in `%Path%` because GCC does not support spaces in filenames. Please follow the next instructions:
	1. Go to Control Panel > System > Advanced system settings > Advanced > Environment Variables.
	2. Find the following environment variables `MSMPI_INC` and `MSMPI_LIB64`.
	3. Open a command window and navigate to `MSMPI_INC`.
	4. Get the *DOS path* of `MSMPI_INC` by executing the following (credit to <https://stackoverflow.com/questions/4051088/how-to-get-dos-path-instead-of-windows-path>): `for %I in (.) do echo %~sI`
	5. Create an environment variable named `MSMPI_INC_GCC` (in Control Panel > System > Advanced system settings > Advanced > Environment Variables). Assign the *DOS path* of `MSMPI_INC` to `MSMPI_INC_GCC`.
	6. Repeat 3-4 for `MSMPI_LIB64`.
	7. Create an environment variable named `MSMPI_LIB64_GCC` (in Control Panel > System > Advanced system settings > Advanced > Environment Variables). Assign the *DOS path* of `MSMPI_LIB64` to `MSMPI_LIB64_GCC`.
- **FICO Xpress 7.8** or superior: Xpress is a commercial solver for mathematical programs. You will need an Xpress license before you can install or download the solver. If you are an academic user, you can request an academic license (free of charge) at <http://subscribe.fico.com/Academic-Partner-Program>.

## Compilation

If you followed the steps above, simply open a Terminal or command line window and execute the compilation script corresponding to your system. If the steps above do not fit your system please contact the [authors](README.md#authors).
