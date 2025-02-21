#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="mergevcf.%J.out"
#SBATCH --job-name="mergevcf"
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

module purge

ml vcftools bcftools


bcftools merge --threads 5 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/02_AD2_n9/AD2.Ah_$seq.n9.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_GD_n25/AD1_GD_n25.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_PR_n43/AD1_PR_n43.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_Yuan_n40/AD1_Yuan_n40.Ah_$seq.combined.vcf.gz \
-Oz -o AD1AD2AD4_n145.Ah_$seq.combined.vcf.gz

bcftools merge --threads 5 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/02_AD2_n9/AD2.Dh_$seq.n9.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_GD_n25/AD1_GD_n25.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_PR_n43/AD1_PR_n43.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_Yuan_n40/AD1_Yuan_n40.Dh_$seq.combined.vcf.gz \
-Oz -o AD1AD2AD4_n145.Dh_$seq.combined.vcf.gz

bcftools view AD1AD2AD4_n145.Ah_$seq.combined.vcf.gz --threads 5 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o AD1AD2AD4_n145.Ah_$seq.combined.bi.vcf.gz

bcftools view AD1AD2AD4_n145.Dh_$seq.combined.vcf.gz --threads 5 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o AD1AD2AD4_n145.Dh_$seq.combined.bi.vcf.gz

module load parallel/20220522-sxcww47
parallel tabix {} ::: AD1AD2AD4_n145.*h_$seq.combined.bi.vcf.gz

