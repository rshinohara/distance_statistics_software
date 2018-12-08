#!/bin/bash
#BSUB -u taki.shinohara@gmail.com      # job name
#BSUB -J t1e_function     # job name
#BSUB -n 50                   # number of tasks in job
#BSUB -R "span[hosts=1]" # run on only one host
#BSUB -q taki_normal              # queue
#BSUB -e "/home/rshi/Distance_Statistics/output/errors.%J.t1e_function"
#BSUB -o "/home/rshi/Distance_Statistics/output/output.%J.t1e_function"

Rscript asympt_sims_Fstat_function_randomG.R
