# Supplemental Materials and Experimental Results for "FPBH.jl: A Feasibility Pump based Heuristic for Multi-objective Mixed Integer Linear Programming in Julia" #

## Understanding the source code: ##

This whole project has been implemented in [Julia v0.6.0]() and divided into six different subpackages: [Modof.jl](https://github.com/aritrasep/Modof.jl), [Modolib.jl](https://github.com/aritrasep/Modolib.jl), [FPBH.jl](https://github.com/aritrasep/FPBH.jl), [CPLEXExtensions.jl](https://github.com/aritrasep/CPLEXExtensions.jl), [FPBHCPLEX.jl](https://github.com/aritrasep/FPBHCPLEX.jl) and [Modoplots.jl](https://github.com/aritrasep/Modoplots.jl) for ease of understanding and maintainence. Let us now describe the purpose and major components of each of these subpackages:

### Modof.jl ###

[Modof.jl](https://github.com/aritrasep/Modof.jl) is a framework used for solving multiobjective mixed integer programs in Julia. Detailed documentation of [Modof.jl](https://github.com/aritrasep/Modof.jl) is available [here](https://aritrasep.github.io/Modof.jl/docs/build/)

1. [ModoModel.jl](https://github.com/aritrasep/Modof.jl/blob/master/src/ModoModel.jl) extends JuMP.jl model for multiple objectives.
2. [Types.jl](https://github.com/aritrasep/Modof.jl/blob/master/src/Types.jl) has efficient data structures for storing instances and solutions of different classes of multiobjective optimization problems.
3. [Utilities.jl](https://github.com/aritrasep/Modof.jl/blob/master/src/Utilities.jl) contains functions for:
    1. selecting, sorting, writing and normalizing a nondominated frontier
    2. computing ideal and nadir points of a nondominated frontier
    3. computing the closest and the farthest point from the ideal and the nadir points respectively.
4. [Quality_of_a_Frontier.jl](https://github.com/aritrasep/Modof.jl/blob/master/src/Quality_of_a_Frontier.jl) contains functions for computing the quality of a nondominated frontier:
    1. exact (for biobjective and triobjective) and approximate (for more than 4 objectives) hypervolumes
    2. cardinality
    3. maximum and average coverage
    4. uniformity
5. [MDLS.jl](https://github.com/aritrasep/Modof.jl/blob/master/src/MDLS.jl) wraps the [MDLS](http://prolog.univie.ac.at/research/MDLS/mdls_code.tar.gz) algorithm for solving multidimensional knapsack and biobjective set packing problems. [MDLS](http://prolog.univie.ac.at/research/MDLS/mdls_code.tar.gz) must be compiled and the respective path of the binaries must be exported as `export PATH="path to mdls binaries:$PATH"`.

### Modolib.jl ###

[Modolib.jl](https://github.com/aritrasep/Modolib.jl) is a collection of instances and their efficient frontiers (if available) of various classes of multiobjective mixed integer programs. It also has function for generating random several classes of random instances. Detailed documentation of [Modolib.jl](https://github.com/aritrasep/Modolib.jl) is available [here](https://aritrasep.github.io/Modolib.jl/docs/build/)

1. [Api.jl](https://github.com/aritrasep/Modolib.jl/blob/master/src/Api.jl) wraps various types of instances and their efficient frontiers (if available) of various classes of multiobjective pure and mixed integer programs.
2. [Generating_Instances.jl](https://github.com/aritrasep/Modolib.jl/blob/master/src/Generating_Instances.jl) contains functions for generating several classes of random instances.

### FPBH.jl ###

[FPBH.jl](https://github.com/aritrasep/FPBH.jl) is the source code of the Feasibility Pump based Heuristic for Multi-objective Mixed Integer Linear Programming. It is developed using [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) and hence supports any LP solver supported by [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl).

1. [starting_solution_creators.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/starting_solution_creators.jl) is the source code of the different weighted sum methods.
2. [feasibility_pumping.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/feasibility_pumping.jl) is the source code of the different feasibility pump methods.
3. [local_search_operators.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/local_search_operators.jl) is the source code of the different local search operators.
4. [decomposition_heuristics.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/decomposition_heuristics.jl) is the source code of stage 1 of FPBH.jl including parallelization.
5. [solution_polishing.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/solution_polishing.jl) is the source code of stage 2 of FPBH.jl including parallelization.
6. In [Overall_Algorithm.jl](https://github.com/aritrasep/FPBH.jl/blob/master/src/Overall_Algorithm.jl) all the different components of FPBH.jl are assembled together.

### CPLEXExtensions.jl ###

[CPLEXExtensions.jl](https://github.com/aritrasep/CPLEXExtensions.jl/blob/master/src/CPLEXExtensions.jl) extends [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) for single-objective optimization by adding additional functionality like deleting constraints, changing coefficient on lhs and rhs of constraints, etc (using `del_constrs!`, `chg_coeffs!`, `set_rhs!`, `chg_coeff_of_obj!`, `chg_coeff_of_rhs!`, `chg_coeff!`, `get_rhs_coef`) and for multi-objective optimization (using `cplex_model`) for [Modof.jl](https://github.com/aritrasep/Modof.jl)

### FPBHCPLEX.jl ###

[FPBHCPLEX.jl](https://github.com/aritrasep/FPBHCPLEX.jl) is the source code of FPBH using [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) and [CPLEXExtensions.jl](https://github.com/aritrasep/CPLEXExtensions.jl). Thus, it only uses CPLEX to solve the underlying LP subproblems. It is important to note that some functions (for generating queues in local search and decomposition heuristics) in [FPBH.jl](https://github.com/aritrasep/FPBH.jl) are reused in [FPBHCPLEX.jl](https://github.com/aritrasep/FPBHCPLEX.jl)

1. [starting_solution_creators.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/starting_solution_creators.jl) is the source code of the different weighted sum methods.
2. [feasibility_pumping.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/feasibility_pumping.jl) is the source code of the different feasibility pump methods.
3. [local_search_operators.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/local_search_operators.jl) is the source code of the different local search operators.
4. [decomposition_heuristics.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/decomposition_heuristics.jl) is the source code of stage 1 of FPBHCPLEX.jl including parallelization.
5. [solution_polishing.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/solution_polishing.jl) is the source code of stage 2 of FPBHCPLEX.jl including parallelization.
6. In [Overall_Algorithm.jl](https://github.com/aritrasep/FPBHCPLEX.jl/blob/master/src/Overall_Algorithm.jl) all the different components of FPBHCPLEX.jl are assembled together.

### Modoplots.jl ###

[Plotting_Nondominated_Frontiers.jl](https://github.com/aritrasep/Modoplots.jl/blob/master/src/Plotting_Nondominated_Frontiers.jl) uses [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) and [Seaborn](https://seaborn.pydata.org/) for plotting nondominated frontiers of 

1. Biobjective:
    1. discrete problems
    2. mixed problems
2. Triobjective problems

## Experiment ##

[MOO_JA_2_Sup.jl](https://github.com/aritrasep/MOO_JA_2_Sup.jl) contains the script for running the whole experiment and generating all the plots. Further, it also contains the [Experimental Results](https://github.com/aritrasep/MOO_JA_2_Sup.jl/blob/master/results/Experimental_Results.csv) and the resulting nondominated frontiers. This repository can be cloned as `git clone https://github.com/aritrasep/MOO_JA_2_Sup.jl`.

### Dependencies: ###

1. [Julia v0.6.0](https://julialang.org/downloads/)
2. [SCIP - v4.0.0](https://github.com/SCIP-Interfaces/SCIP.jl)
3. [Gurobi - v7.5](https://github.com/JuliaOpt/Gurobi.jl)
4. [CPLEX - v12.7](https://github.com/JuliaOpt/CPLEX.jl)
5. [FPBHCPLEX.jl](https://github.com/aritrasep/FPBHCPLEX.jl)
6. [MDLS](http://prolog.univie.ac.at/research/MDLS/mdls_code.tar.gz) must also be compiled and the respective path of its binaries must be exported as `export PATH="path to mdls binaries:$PATH"`.

### Running the whole Experiment: ###

1. Install `Julia v0.6.0`, `SCIP`, `Gurobi` `CPLEX`, `FPBHCPLEX.jl` and `MDLS`
2. Start a `julia` terminal inside the `src` directory with atleast 4 threads using `julia -p 4`
3. Type the following in the terminal to execute the whole experiment

```julia
include("Experimental_Run.jl")
run_experiment()
```

The above commands using [Experimental_Run.jl](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/src/Experimental_Run.jl), will run all experiments related to MDLS, V1, V2, V3, V4 and V5 T1 in parallel using 4 ( or more ) workers. However, all experiments related to V5 T2, V5 T3 and V5 T4 will be run serially.

### Experimental Results: ###

All results will be written as csv files inside `results` directory. [Results](https://github.com/aritrasep/MOO_JA_2_Sup.jl/blob/master/results/Experimental_Results.csv) is the detailed experimental results for our paper.

### Nondominated Frontiers: ###

1. Nondominated Frontiers:
	1. [V1](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V1/)
	2. [V2](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V2/)
	3. [V3](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V3/)
	4. [V4](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V4/)
	5. V5:
		1. [1 Thread](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V5%20T1/)
		1. [2 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V5%20T2/)
		1. [3 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V5%20T3/)
		1. [4 Threads](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/V5%20T4/)
	6. [True / Reference Frontiers](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/results/non_dom_sols/True_Frontier/)

### Summarizing Results and Generating Plots: ###

**Modoplots.jl** must be installed on the local machine. If it has not been installed, it can be done so using the following instructions in a **Julia** terminal:

```julia
Pkg.clone("https://github.com/aritrasep/Modoplots.jl")
Pkg.build("Modoplots")
```

Start a `Julia` terminal inside the `src` directory and type the following in the terminal to generate plots used in the article

```julia
include("Summarizing_Results.jl")
```

The above commands using [Summarizing_Results.jl](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/src/Summarizing_Results.jl), will generate [all plots](https://github.com/aritrasep/MOO_JA_2_Sup.jl/tree/master/plots/) inside the `plots` directory.
