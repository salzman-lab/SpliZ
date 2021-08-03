#!/bin/bash
#
#SBATCH --job-name=nextflow
#SBATCH --output=job_output/nextflow.%j.out
#SBATCH --error=job_output/nextflow.%j.err
#SBATCH --time=24:00:00
#SBATCH -p horence,quake
#SBATCH --nodes=1
#SBATCH --mem=60Gb

conf=$1
nextflow run main.nf -c ${conf} -resume
