[![DOI](https://zenodo.org/badge/245032003.svg)](https://zenodo.org/badge/latestdoi/245032003)

# AsyncLSD: Asynchronous Lagrangian Scenario Decomposition

Implementents asynchronous Lagrangian scenario decomposition for two-stage mixed integer programs with (weak) relatively complete recourse. 

## What's this?

This repository contains the source code of the implementation of the asynchronous decomposition algorithm for stochastic mixed-integer programs developed/presented in [this research article](#citation). The algorithm implementation takes a two-stage stochastic program with (weak) relatively complete recourse in SMPS format as input. It computes a feasible solution and an optimimality gap associated with that solution. The implementation uses FICO Xpress to solve subproblems, called through its C API, and MPI to manage communications between subprocesses.

## Usage

Please follow the instructions in [INSTALLATION.md](INSTALLATION.md) in order to compile the source code into an executable for your system.

To execute the program in parallel in interactive mode (e.g. in your local machine or at an assigned node of a cluster) simply use

```
mpiexec -n $num_procesors AsyncLSD $WorkDir $InstanceDir $Instance $OptionsFile
```

where:

- `$num_procesors`:  Number of MPI processes.
- `$WorkDir`:        Working directoy. Here we will write all files concerning execution of the algorithm.
- `$InstanceDir`:    Directory containing the instance in SMPS, compressed as `$Instance.tar.gz` (containing `$Instance.cor`, `$Instance.sto` and `$Instance.tim`)
- `$Instance`:       Name of the instance to be solved.
- `$OptionsFile`:    File describing the configuration parameters for the execution of the decomposition algorithm and for the solution of subproblems using Xpress. See [PARAMETERS.md](PARAMETERS.md) for further details on how to create an option file.

In High Performance Computing (HPC) clusters you might need to replace `mpiexec` by a command specific to the job scheduler managing the HPC cluster. Please consult the documentation of the job scheduler or contact the [authors](#authors) if you require help.

## Authors

* [Ignacio Aravena](https://sites.google.com/site/iaravenasolis/) -- algorithm design and implementation
* [Anthony Papavasiliou](https://perso.uclouvain.be/anthony.papavasiliou/public_html/) -- algorithm design

## Citation

If you find AsyncLSD useful in your work, we kindly request that you cite the following paper:

```
@article{AravenaPapavasiliou2020,
author = {Ignacio Aravena and Anthony Papavasiliou},
title = {Asynchronous Lagrangian Scenario Decomposition},
journal = {Mathematical Programming Computation},
volume = {},
number = {},
pages = {},
year = {2020},
doi = {},
}
```

## Acknowledgement

The development of this algorithm was partially funded by the ENGIE Chair on Energy Economics and Energy Risk Management, by the Univeriste catholique de Louvain through an FSR grant, and by the U.S. Department of Energy through the Lawrence Livermore National Laboratory under contract DE-AC52- 07NA27344.
