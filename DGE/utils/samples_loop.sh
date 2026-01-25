#!/bin/bash

# Parameters
NAME="LINENAME"
max_index=27  # N. files starting as NAME
stranded="reverse"

# Sample loop
for ((i=1; i<=$max_index; i++)); do

    sbatch submit_slurm_job.sh $NAME $i $stranded

done
