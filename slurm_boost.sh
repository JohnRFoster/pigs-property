#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=8               # for a multithredded job
#SBATCH --exclusive
#SBATCH --mem=48g                       # memory
#SBATCH --job-name=oosBoost                 # job name
#SBATCH --output=outfiles/oosBoostMean.txt      # output file

module add R
Rscript 2_workflow_ml.R $SLURM_ARRAY_TASK_ID
