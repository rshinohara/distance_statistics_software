#!/bin/bash
#BSUB -u taki.shinohara@gmail.com      # job name
#BSUB -J cure_sims[1-48]      # job name
#BSUB -n 15                  # number of tasks in job
#BSUB -R "span[hosts=1]" # run on only one host
#BSUB -q taki_normal             # queue
#BSUB -e "/home/rshi/Distance_Statistics/output/errors.%J.cure_sims"
#BSUB -o "/home/rshi/Distance_Statistics/output/output.%J.cure_sims"

Rscript cure_simulations.R
