#!/bin/bash
#BSUB -u taki.shinohara@gmail.com      # job name
#BSUB -J timeapprox_power_graph      # job name
#BSUB -n 60                   # number of tasks in job
#BSUB -R "span[hosts=1]" # run on only one host
#BSUB -q taki_normal              # queue
#BSUB -e "/home/rshi/Distance_Statistics/output/errors.%J.timeapprox_power_graph"
#BSUB -o "/home/rshi/Distance_Statistics/output/output.%J.timeapprox_power_graph"

Rscript asympt_sims_Fstat_graph_randomG_power_findapprox.R