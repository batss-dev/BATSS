# Parallel and cluster computation

To speed up computation, the `batss.glm` function can use
parallelisation on a single machine by setting
`computation = "parallel"` (the default; see the help of `batss.glm` for
details).

When using multiple machines and/or a cluster, parallelisation can be
achieved by splitting the set of seeds - corresponding to as many
simulated trials - between machines or cluster CPUs, saving the
`batss.glm` outputs, and merging them using the function
`batss.combine`.

We describe here how to use the functions `batss.glm` and
`batss.combine` in two settings when:

- using several machines,
- using a cluster.

## Several machines

Let’s assume a **BATSS** user wants to perform a Monte Carlo simulation
considering 10,000 trials and has two computers, each of them with 10
CPUs.

The strategy we suggest consists of

- running `batss.glm` on each each computer for a different subset of
  5,000 seeds among the 10,000 seeds of interest specified in the
  argument `R` (so that each computer evaluates a different set of
  seeds), with
  - argument `computation` set to `parallel` and the argument `mc.cores`
    set to `10` (or
    [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)),
    the number of CPUs of the computer of interest,
  - argument `extended` set to `1` or `2` (check
    [`?batss.glm`](../../reference/batss.glm.md) for details),
- saving the `batss.glm` output as RData files with the function `save`
  under a name specific to the task of each computer (such as the first
  or last seed evaluated by each computer),
- finally, using the function `batss.combine` to combine these outputs.

Let’s have the first computer run seeds 1 to 5,000 and save the output
under `seed1to5000.rdata`:

``` r
require(BATSS)

##
## run seeds 1 to 5000 (example of batss.glm)
##

seed1to5000 = batss.glm(
                model            = y ~ group,   
                var              = list(y     = rnorm,
                                        group = alloc.balanced),
                var.control      = list(y = list(sd = 5)),
                beta             = c(1, 1, 2),
                which            = c(2:3),
                alternative      = "greater",
                R                = 1:5000,
                N                = 200,
                interim          = list(recruited = seq(100, 180, 20)),
                prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
                eff.arm          = eff.arm.simple,
                eff.arm.control  = list(b = 0.975),
                fut.arm          = fut.arm.simple,
                fut.arm.control  = list(b = 0.05),
                computation      = "parallel",
                H0               = TRUE,
                mc.cores         = 10,
                extended         = 1)

##
## save results as an rdata file
##

# specify here the folder in which to save the output
path_computer1 = "~/"
# save
save(seed1to5000, file = paste0(path_computer1,"seed1to5000.rdata"))
```

Let’s now have the second computer run seeds 5,001 to 10,000 and save
the output under `seed1to5000.rdata`:

``` r
##
## run seeds 5001 to 10000 (example of batss.glm)
##

seed5001to10000 = batss.glm(
                model            = y ~ group,   
                var              = list(y     = rnorm,
                                        group = alloc.balanced),
                var.control      = list(y = list(sd = 5)),
                beta             = c(1, 1, 2),
                which            = c(2:3),
                alternative      = "greater",
                R                = 5001:10000,
                N                = 200,
                interim          = list(recruited = seq(100, 180, 20)),
                prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
                eff.arm          = eff.arm.simple,
                eff.arm.control  = list(b = 0.975),
                fut.arm          = fut.arm.simple,
                fut.arm.control  = list(b = 0.05),
                computation      = "parallel",
                H0               = TRUE,
                mc.cores         = 10,
                extended         = 1)

##
## save results as an rdata file
##

# specify here the folder in which to save the output
path_computer2 = "~/"
# save
save(seed5001to10000, file = paste0(path_computer2,"seed5001to10000.rdata"))
```

Transfer the objects `seed1to5000.rdata` and `seed5001to10000.rdata` to
the same folder of the same computer, for example `path_computer1`
above, and merge the objects with the function `batss.combine` as
follows:

``` r
# combine
seed1to10000 = batss.combine(
    paths = paste0(path_computer1, c("seed1to5000","seed5001to10000"),".rdata"))

# look at combined results
summary(seed1to10000)
```

## Cluster

Let’s assume a **BATSS** user wants to perform a Monte Carlo simulation
considering 10,000 trials and has access to 500 CPUs of a cluster.

The strategy we suggest consists of

- running `batss.glm` on each CPU for a specific subset of the 10,000
  seeds of interest specified in the argument `R` (so that each CPU
  evaluates a different set of seeds), with argument `computation` set
  to `sequential` (as parallelisation is already achieved by the large
  number of CPUs),
- saving each CPU’s `batss.glm` output as an RData file using the
  function `save` with filenames indicating the task (e.g., based on the
  first or last seed evaluated by each CPU),
- using the `batss.combine` function to merge all these outputs into a
  single one.

There are multiple ways to accomplish this. In the following, we
describe the approach we follow. Let us assume that

- the cluster runs a **Linux OS** with **Slurm** as workload manager (a
  common setup in cluster computing),
- the working directory for this simulation is `~/batss-example/`, which
  contains:
  - a folder `~/batss-example/in/` with (optional) inputs, like
    parameters values saved in RData format, or an R script containing
    user-defined functions to be used by `batss.glm` (like the function
    `treatalloc.fun` described in the ANCOVA help page) that is to be
    sourced by each CPU,
  - a folder `~/batss-example/out/` where all outputs from the CPUs will
    be stored,
  - an R script `~/batss-example/1-run.r`, which
    - optionally loads and sources elements of the folder
      `~/batss-example/in/`,
    - runs the `batss.glm` function for a set of seeds assigned by
      **Slurm**,
    - saves the output in folder `~/batss-example/out/` folder,
  - an R script `~/batss-example/2-combine.r`, which merges the
    different 500 `batss.glm` outputs with the function `batss.combine`,
  - a shell script `~/batss-example/batss_sim.sh` containing the
    instructions for each CPU.

The simulation is conducted by:

- i/ invoking the `sbatch` command, which
- ii/ runs the shell script `~/batss-example/batss_sim.sh` on each CPU,
  which
- iii/ starts R and executes the `~/batss-example/1-run.r` script that
  performs the parallelised simulation,
- iv/ and finally executing the `~/batss-example/2-combine.r` script
  that merges all simulation results.

Let’s describe each step:

### sbatch command

The following command corresponds to the **Slurm** `sbatch` submission
command. It submits a job array with task IDs from 1 to 500 (i.e., the
number of CPUs), executing the script sim_batch.sh and passing `in`
(input folder) and `out` (output folder) as arguments for each task:

``` bash
# move to the folder of interest 
cd ~/batss-example
# sbatch function
sbatch --array=1-500 batss_sim.sh in out
```

For additional options such as setting the maximum computation time,
memory limits, or email notifications, please refer to the help section
of the sbatch command or consult your cluster’s documentation.

### batss_sim.sh shell script

The following command of the shell script fed into the `sbatch` call
above tells each CPU to start `R` (both in `slave` and `vanilla` mode),
run the script `1-run.r` with

- arguments `$1` and `$2`, corresponding respectively to the `in` and
  `out` folders specified at the end of the call to `sbatch`,
- the location where to save the `Rout` files related to each task
  (i.e., CPU job here referred to as`$SLURM_ARRAY_TASK_ID` and provided
  by **Slurm**): these files will be saved in the `out` folder
  (indicated as `$2`).

``` bash
R --slave --vanilla < 1-run.r --args $SLURM_ARRAY_TASK_ID $1 $2 > $2/${SLURM_ARRAY_TASK_ID}.Rout 2>&1
```

Note that the R program may not be directly available via a call to R
and that you might need to

- specify the full path to R,
- add to the shell script a way to make R available (like the `module`
  command, for example)

You can check this with your cluster manager.

### 1-run.r R script

The following code shows the content of the R script run by each CPU.
The code

- loads the library **BATSS**,
- defines the vector of seeds related to the job ID attributed by
  **Slurm**,
- runs the batss.glm for that vector of seeds with `computation` set to
  `sequential` and `extended` set to 1,
- saves the results in the output folder under of name that corresponds
  to the first seed.

``` r
###################################
## setup 
###################################


# load library
library(BATSS)

# define list of seed for CPU to handle where
# - args[1] is the job id attributed by slurm 
#   that we will use as first seed
# - 10000 is the number of simulations (i.e., seeds)
# - 500 is the number of CPUs
seed.list = seq(as.numeric(args[1]),10000,500)

# define path to in/ and out/ folders (input of sbatch)
path_in  <- paste0(args[2],"/")
path_out <- paste0(args[3],"/")
# optionally load/source info from relevant files 
# of "in/" here


###################################
## simulation 
###################################

# only compute res if needed
if(!any(dir(path_out)==id.seed$id[1])){
    ##
    ## batss.glm
    ##
    start = Sys.time()                        
    sim   = batss.glm(
            model            = y ~ group,   
            var              = list(y     = rnorm,
                                    group = alloc.balanced),
            var.control      = list(y = list(sd = 5)),
            beta             = c(1, 1, 2),
            which            = c(2:3),
            alternative      = "greater",
            R                = seed.list,
            N                = 200,
            interim          = list(recruited = seq(100, 180, 20)),
            prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
            eff.arm          = eff.arm.simple,
            eff.arm.control  = list(b = 0.975),
            fut.arm          = fut.arm.simple,
            fut.arm.control  = list(b = 0.05),
            computation      = "sequential",
            H0               = TRUE,
            extended         = 1)
    finish = Sys.time()

    ## print required time
    cat("\trequired time: ",finish-start,"\n\n")

    ## store results
    save(sim,file=paste0(path_out,seed.list[1],".rdata"))
}# end if

cat("\n\t",date(),"\n")
cat("\n\t DONE!\n")
q("no")
```

### 2-combine.r R script

Once the simulation is complete, the following code merges all results
into a single object:

``` r
# load library
library(BATSS)

# define list of 'successful' jobs
job.list = dir("out/")[!grepl("Rout",dir("out/"))]

# list of potential 'unsuccessful' slurm jobs to be 
# investigated by looking at the corresponding 
# 'out' (slurm) and 'Rout' (R) files
seq(1,500)[is.na(match(paste0(seq(1,500),".rdata"),job.list))]

# merge
sim = batss.combine(paste0("out/",job.list))

# store results
save(sim,file=paste0("out/full.rdata"))
```
