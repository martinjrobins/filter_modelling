#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

#SBATCH --cpus-per-task=1

# set max wallclock time
#SBATCH --time=100:00:00

# set name of job
#SBATCH --job-name=filter

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=martin.robinson@cs.ox.ac.uk

# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set

module load python/anaconda2/4.3.0
source activate rnb
source job_vis_${SLURM_ARRAY_TASK_ID}.sh
~
~
