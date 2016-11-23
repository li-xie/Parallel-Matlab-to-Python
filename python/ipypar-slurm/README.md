# Setting up ipyparallel @hutch (Slurm Version)

## Setup the Environment

There seems to be some inconsistencies in the configuration of the ipyparallel
module in `/app`- we should also be using the most recent versions of the tools
(ipython et alia) as ipyparallel is going through a refactor.

Thus, at the moment, it is necessary to configure a few things manually before
starting the ipyparallel cluster

 - `module load python3/3.5.0`
 - `pip3 install --user ipython3`
 - `pip3 install --user ipyparallel`
 - `export PATH=~/.local/bin:${PATH}`

## Start the Cluster

 - copy the `profile_slurm` directory in this directory to `~/.ipython/`.  This
   will make the profile available via the `--profile=` option
 - start the cluster using: 
   `ipcluster start -n ${SLURM_NTASKS} --profile='slurm'`

In this instance, ipyparallel takes over and submits the job to slurm.  You'll
see a message indicating the job ID.  The cluster will be ready when you see
the message:

    2016-11-22 15:29:20.863 [IPClusterStart] Engines appear to have started successfully

Note that the ipcluster

You can then run the sample script: `./run_sim.py`.  If successful you will get
a mean of the `c` values

