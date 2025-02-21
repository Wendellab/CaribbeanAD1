#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G 
#SBATCH --time=05:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.plink.%J.out"
#SBATCH --job-name="plink"

ml bcftools

module load plink/1.90b6.21

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/03_AD1AD2AD4_n145/AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf
output=AD1AD2AD4_n145

#Filtering LD by Plink 
plink --threads 30 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 30 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

plink --threads 30 --vcf $vcf --allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
