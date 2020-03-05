# Algorithm parameters

To create a customized configuration parameter file use the included example file [example_options_file.txt](example/example_options_file.txt) as a base and modify it. Each mandatory field and its possible values are described below:

| Parameter | Values | Description |
| :---      | :----: | :---        |
| `VERBOSITY` | 0,..,4 | Verbosity level, in increasing order |
| `MAXITER` | 1,... | Maximum number of iterations normalized by the number of scenarios, i.e. `MAXITER`=100 sets a maximum of 100*N* coordinate descent iterations |
| `MAXTIME` | (0, +Inf) | Maximum execution time in seconds |
| `DSHARE_MIDITER`<sup>1</sup> | (0, 1) | Medium point for dynamic change in resource allocation policy for processors in terms of the total number of iterations |
| `DSHARE_START`<sup>1</sup> | (0, 1) | Proportion of processors allocated to solve dual subproblems in the coordinate-descent first iteration |
| `DSHARE_END`<sup>1</sup> | (0, 1) | Proportion of processors allocated to solve dual subproblems at the coordinate-descent last iteration, `MAXITER`*N* |
| `NONSCEN_PERIOD` | (0, 1) | Periodicity to solve subproblem (6), e.g. `NONSCEN_PERIOD`=0.5 would cause that subproblem (6) is solved every two subproblems (7) |
| `NONSCEN_CENTER`<sup>2</sup> | 10, 11 | See *u*<sub>0</sub> in eq. (8). 10: fixed *u*<sub>0</sub> at the average of the first-stage solutions during the first evaluation for each scenario, 11: use incumbent first stage solution as *u*<sub>0</sub> |
| `DUAL_STEPSIZE` | 20, 21, 22 | 20: constant, 21: diminishing, 22: Polyak |
| `DUAL_STEPSIZE_P` | (0, +Inf) | *p*<sub>0</sub> in eq. (14) |
| `DUAL_STEPSIZE_R` | (0, +Inf) | *r* in eq. (14) |
| `DUAL_STEPSIZE_Q` | (0.5, 1] | *q* in eq. (14) |
| `DUAL_STEPSIZE_XI` | (0, 1) | &xi; in eq. (14), expressed in relative terms with respect to *LB*<sub>*k*</sub> |
| `DUAL_STEPSIZE_SGM` | (0, +Inf) | &sigma; in eq. (14) |
| `DUAL_GNORM_FILTER` | (0, 1) | Allows to put a filter in the updating of *g*<sub>*k*</sub> as new information arrives from subproblems. `DUAL_GNORM_FILTER`=0 implies no filter, `DUAL_GNORM_FILTER`=1 implies that *g*<sub>*k*</sub> is not updated. |
| `PRIMAL_RECOVERY` | 30,..,33 | 30: FIFO, 31: Random, 32: LIFO, 33: IS |
| `PRIMAL_IS_SSIZE` | (0, 1) | Relative size of subsample to form average when using IS primal recovery |
| `INIT_TYPE` | 40, 41 | 40: run initialization using LP relaxations, 41: run initialization using period relaxations |
| `DUAL_GAP_TOL` | (0, 1) | Duality gap tolerance for termination; `DUAL_GAP_TOL`=0.01 would terminate the algorithm once the duality gap reaches 1% |
| `SMOOTHING_PARAM`<sup>2</sup> | (0, +Inf) | &mu; in eq. (8) |
| `PROJ_ROUND_MIDP` | (0, 1) | Conservativeness of rounding operations. Variables with non-integer variables will tend to be rounded to 1 if they are `PROJ_ROUND_MIDP` or more. Used during initialization with period relaxation and by IS primal recovery. |
| `PERIOD_AGG_POL` | 50, 51, 52 | Policy for aggregating values of variables appearing in constraints in multiple periods when performing initialization with period relaxation. 50: use minimum value, 51: use average, 52: use maximum value |
| `SCEN_LP_ALG` | p,d,b | Algorithm to solve the root LP relaxation of the scenario subproblem. p: primal simplex, d: dual simplex, b: barrier |
| `PERIOD_LP_ALG` | p,d,b | Algorithm to solve the root LP relaxation of the period subproblem (initialization). p: primal simplex, d: dual simplex, b: barrier |

<sup>1</sup> The proportion of the processors dedicated to dual iterations follows an hyperbolic tangent with *DualShare*=`DSHARE_START` at iteration 1, *DualShare=*`DSHARE_MIDITER` at iteration `MAXITER`*N*/2, and *DualShare*=`DSHARE_END` at iteration `MAXITER`*N*.

<sup>2</sup> For performance of the method and accuracy of the dual function smooth approximation, we recommend setting `NONSCEN_CENTER`=11 and `SMOOTHING_PARAM` to a small value for your problem, e.g. `SMOOTHING_PARAM`=0.1.

The implementation of AsyncLSD, additionally, allows to pass any configuration parameter to Xpress using their code numbers (see `Xpress.h`). These parameters are optional and are declared in the following sections of the configuration file:
	
| Section | Parameters |
| :---    | :---       |
| `XPRESS_SCEN_INT` | Integer parameters for scenario subproblems |
| `XPRESS_SCEN_DBL` | Double parameters for scenario subproblems |
| `XPRESS_PERIOD_INT` | Integer parameters for period-decomposition subproblems and first-stage subproblem |
| `XPRESS_PERIOD_DBL` | Double parameters for period decomposition subproblems and first stage subproblem |
