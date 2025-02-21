#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G 
#SBATCH --time=15:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="plink.%J.out"
#SBATCH --job-name="plink"


ml vcftools bcftools

output=MKGDPR_n86

bcftools view -S subset_n86.txt /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/03_AD1AD2AD4_n145/AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor \
-o $output.AhDh.combined.vcf


##########################################

module load plink/1.90b6.21

vcf=$output.AhDh.combined.vcf

#Filtering LD by Plink 
plink --threads 10 --vcf $vcf  \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 10 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

plink --threads 10 --vcf $vcf --allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --genome --pca 20 var-wts --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R


module purge

module load py-seaborn/0.12.2-py310-6p2ciw6

python 3D-PCA-plot-wx3.py --evec GD_n90.eigenvec.3d --eval GD_n90.eigenval.3d --s 100 --x 8 --y 8 --o MKGDPR_n86_3dPCA.tiff
