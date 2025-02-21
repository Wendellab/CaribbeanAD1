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
-O MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename.MKAD2AD4AD5_n51.txt MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz -o MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf >  MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf

grep -v "#" MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf | wc -l

