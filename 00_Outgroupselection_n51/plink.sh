#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G 
#SBATCH --time=05:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.plink.%J.out"
#SBATCH --job-name="plink"

ml bcftools

module load plink/1.90b6.21

#Filtering LD by Plink 
plink --threads 30 --vcf ../MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf  \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out MKAD2AD5AD4_n51

plink --threads 30 --extract MKAD2AD5AD4_n51.prune.in  \
--make-bed --allow-extra-chr --const-fid --out MKAD2AD5AD4_n51 \
--recode vcf-iid --vcf ../MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf

plink --threads 30 --vcf ../MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf --allow-extra-chr --const-fid \
--extract MKAD2AD5AD4_n51.prune.in --recode \
--make-bed --pca 20 var-wts --distance square 1-ibs --out MKAD2AD5AD4_n51

paste -d '\t'  MKAD2AD5AD4_n51.mdist.id MKAD2AD5AD4_n51.mdist >  MKAD2AD5AD4_n51.distmatrix

module load r

Rscript PCA_plot.R
