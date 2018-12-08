#!/bin/bash
#BSUB -u taki.shinohara@gmail.com      # job name
#BSUB -J t1e_graph      # job name
#BSUB -n 50                   # number of tasks in job
#BSUB -R "span[hosts=1]" # run on only one host
#BSUB -q cceb_normal              # queue
#BSUB -e "/home/rshi/Distance_Statistics/output/errors.%J.t1e_graph"
#BSUB -o "/home/rshi/Distance_Statistics/output/output.%J.t1e_graph"

Rscript asympt_sims_Fstat_graph_randomG.R