#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=7               # for a multithredded job
#SBATCH --mem=96g                       # memory
# #SBATCH --time=144:00:00
#SBATCH --job-name=iteratuveBatch             # job name
#SBATCH --output=outfiles/iterativeBatch.txt    # output file

module add R
Rscript 1_workflow_mcmc.R
