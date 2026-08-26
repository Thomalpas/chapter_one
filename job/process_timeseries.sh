#!/bin/bash
#
#SBATCH --job-name=process_timeseries
# Set number of iteration
#SBATCH --array=1-90%90
# Amount of RAM requested per job
#SBATCH --mem=64G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
# Request
#SBATCH --time=2-00:00:00

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/julia_coding/chapter_one
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

JULIA=/users/bi1tma/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

STEPSIZE=1000

START=$((1 + ${STEPSIZE} * (${SLURM_ARRAY_TASK_ID} - 1)))
END=$((${STEPSIZE} * ${SLURM_ARRAY_TASK_ID}))

echo "Starting task from ${START} to ${END}"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/processing_timeseries.jl ${START} ${END} "HPC"
