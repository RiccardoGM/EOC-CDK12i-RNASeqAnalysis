#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH --output=htseq_count_%A.out
#SBATCH --error=htseq_count_%A.err
#SBATCH --partition=regular1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # Adjust as needed
#SBATCH --mem=32G          # Adjust as needed
#SBATCH --time=12:00:00    # Adjust as needed, example (2hr): --time=2:00:00


# --- Print info --- #
#
NOW=`date +%H:%M-%a-%d/%b/%Y`
echo '------------------------------------------------------'
echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
echo 'Job is running on node(s): '
echo  $SLURM_JOB_NODELIST
echo '------------------------------------------------------'
echo 'WORKINFO:'
echo 'SLURM: job starting at           '$NOW
echo 'SLURM: sbatch is running on      '$SLURM_SUBMIT_HOST
echo 'SLURM: executing on cluster      '$SLURM_CLUSTER_NAME
echo 'SLURM: executing on partition    '$SLURM_JOB_PARTITION
echo 'SLURM: working directory is      '$SLURM_SUBMIT_DIR
echo 'SLURM: current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
echo ""
echo 'JOBINFO:'
echo 'SLURM: job identifier is         '$SLURM_JOBID
echo 'SLURM: job name is               '$SLURM_JOB_NAME
echo ""
echo 'NODEINFO:'
echo 'SLURM: number of nodes is        '$SLURM_JOB_NUM_NODES
echo 'SLURM: number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
echo 'SLURM: number of gpus/node is    '$SLURM_GPUS_PER_NODE
echo '------------------------------------------------------'


# --- Activate Conda --- #
. "/path/to/conda.sh"
conda activate htseq_env


# --- Retrieve inputs from command line --- #
NAME="$1"
idx="$2"
stranded="$3" # strandedness


# --- Run the Python script --- #
python count_reads.py "$NAME" "$idx" "$stranded"


# --- Print a job completion message --- #
echo "Job completed!"
