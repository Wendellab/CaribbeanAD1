#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="Picard.out"
#SBATCH --job-name="picard"

module load picard/2.27.4

picard GatherVcfs \
$(for vcf in *.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O AD1AD2AD4_n145.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename.AD1AD2AD4_n145.txt AD1AD2AD4_n145.AhDh.combined.bi.vcf.gz -o AD1AD2AD4_n145.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" AD1AD2AD4_n145.AhDh.combined.bi.rehead.vcf >  AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf

grep -v "#" AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf | wc -l