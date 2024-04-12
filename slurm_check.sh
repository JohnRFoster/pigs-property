#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=2               # for a multithredded job
#SBATCH --mem=24g                       # memory
#SBATCH --job-name=check_300                 # job name
#SBATCH --output=outfiles/check_300.txt      # output file

module add R
Rscript R/check_mcmc.R
