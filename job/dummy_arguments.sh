#!/bin/bash
#
#SBATCH --job-name=install_packages
# Amount of RAM requested per job
#SBATCH --mem=8G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/julia_coding/chapter_one
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

JULIA=/users/bi1tma/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/parse_argument.jl "HPC"