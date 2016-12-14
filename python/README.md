# Setting up ipyparallel @hutch

## Setup the Environment

There seems to be some inconsistencies in the configuration of the ipyparallel
module in `/app`- we should also be using the most recent versions of the tools
(ipython et alia) as ipyparallel is going through a refactor.

Thus, at the moment, it is necessary to configure a few things manually before
starting the ipyparallel cluster

```
module load python3/3.5.0
pip3 install --user ipython3
pip3 install --user ipyparallel
export PATH=~/.local/bin:${PATH}
```

## Starting the Cluster

### Using MPI in a Slurm Allocation

 - copy `ipypar-mpi/ipyparallel/profile_mpi` to `~/.ipython`. This enables the
   MPI engines and puts in some sensible defaults
 - get an allocation from Slurm using `salloc -n <tasks>`
 - once you've received the allocation, start the cluster with
   `ipcluster start -n ${SLURM_NTASKS} --profile='mpi'`

The cluster is ready when you see the message:

    2016-11-22 15:29:20.863 [IPClusterStart] Engines appear to have started successfully

### Using the Slurm Controller and Engine

 - copy `ipypar-slurm/ipyparallel/profile_slurm` to `~/.ipython` to enable
   the slurm engine with sensible defaults
 - start the cluster using: `ipcluster start -n <tasks> --profile='slurm'`

In this instance, ipyparallel takes over and submits the job to slurm.  You'll
see a message indicating the job ID.  The cluster will be ready when you see
the message:

    2016-11-22 16:21:57.268 [IPClusterStart] Job submitted with job id: '45374171'
    2016-11-22 15:29:20.863 [IPClusterStart] Engines appear to have started successfully

## Running Tasks

The first thing you will notice (assuming successful completion of either of the above cluster-startups) is that these controllers run in the foreground.  Running tasks with the `ipyparallel` framework requires that your task can access the controller processes.  There are a couple ways of doing this:

 - run `ipcluster` with `--daemonize` to put it in the background
 - background `ipcluster` with `&` or `^Z` and `bg`
 - start another shell to the _same host_ (e.g. if you start `ipcluster` on
   rhino2, open a new shell to rhino2
 - configure the ipyparallel client to talk to the IP where you started
   `ipcluster` (I don't know how to do this yet)

Once you have a shell, the cluster is available for you to run your tasks.  The
`run_sim` commands in each of the subdirectories is such a command that can be
run on these clusters.
