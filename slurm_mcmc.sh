#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=3               # for a multithredded job
#SBATCH --mem=50g                       # memory
#SBATCH --job-name=dev3                 # job name
#SBATCH --output=outfiles/dev3.txt      # output file

module add R
Rscript R/1_workflow_mcmc.R SLURM_CPUS_PER_TASK
