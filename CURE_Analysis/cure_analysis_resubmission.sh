#!/bin/bash
#BSUB -u taki.shinohara@gmail.com      # job name
#BSUB -J cure_analysis      # job name
#BSUB -n 1                  # number of tasks in job
#BSUB -R "span[hosts=1]" # run on only one host
#BSUB -q taki_normal             # queue
#BSUB -e "/home/rshi/Distance_Statistics/output/errors.%J.cure_analysis"
#BSUB -o "/home/rshi/Distance_Statistics/output/output.%J.cure_analysis"

Rscript cure_analysis_resubmission.R
