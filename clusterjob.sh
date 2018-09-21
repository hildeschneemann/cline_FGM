#!/bin/bash

# Shell to use
#$ -S /bin/bash
# Name of the job in SGE
#$ -N sim_original_cline
# Name of the queue to use
#$ -q small.q
# Maximum hardware time allowed for this job
#$ -l h_rt=10:00:00
# run in the current directory
#$ -cwd

module load singularity-2.5.1

singularity exec ubuntu_Hdf5_Boost.simg ./a.out
