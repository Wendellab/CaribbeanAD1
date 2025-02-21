#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --exclude=speedy[9-10],speedy6
#SBATCH --open-mode=append
#SBATCH --output="mergevcf.%J.out"
#SBATCH --job-name="mergevcf"
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

ml vcftools bcftools

bcftools view -S AD2_n9.txt /work/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Ah_$seq.combined.vcf.gz -O z -o AD2.Ah_$seq.n9.combined.vcf.gz
bcftools view -S AD2_n9.txt /work/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Dh_$seq.combined.vcf.gz -O z -o AD2.Dh_$seq.n9.combined.vcf.gz

module load parallel/20220522-sxcww47
parallel tabix {} ::: AD2.*h_$seq.n9.combined.vcf.gz