#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G 
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n90.%J.out"
#SBATCH --job-name="jointGeno_n90"


module load r

Rscript LEA.R

cut -d ' ' -f 2 *.ped > samplename.txt