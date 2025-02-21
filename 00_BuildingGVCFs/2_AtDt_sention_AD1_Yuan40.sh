#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --exclude=speedy[9-10],speedy6
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_Yuan_n40.%J.out"
#SBATCH --job-name="YuanVCF_n40"
#SBATCH --array=1-13


seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

Dir=/work/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval
ref=/work/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=20 #NUMBER_THREADS

output1=AD1_Yuan_n40

cd $Dir
#mkdir $output1

module load sentieon-genomics/202308-mdtz2zq
export SENTIEON_LICENSE=reimu.las.iastate.edu:8990

#joint SNP calling:

cat AD1_Yuan_n40.txt | sentieon driver --interval Ah_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Ah_$seq.vcf -
mv $TMPDIR/$output1.Ah_$seq.vcf* $Dir/$output1/

cat AD1_Yuan_n40.txt | sentieon driver --interval Dh_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Dh_$seq.vcf -
mv $TMPDIR/$output1.Dh_$seq.vcf* $Dir/$output1/

#echo "Filtering VCF file using vcftools and rename samples using bcftools"
cd $Dir/$output1/

ml vcftools bcftools

vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Ah_$seq.variant
vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Ah_$seq.invariant

vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Dh_$seq.variant
vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Dh_$seq.invariant

module load parallel/20220522-sxcww47

parallel bgzip {} ::: $output1.*h_$seq.*variant.recode.vcf
parallel tabix {} ::: $output1.*h_$seq.*variant.recode.vcf.gz

bcftools concat --allow-overlaps --threads $thr $output1.Ah_$seq.variant.recode.vcf.gz $output1.Ah_$seq.invariant.recode.vcf.gz -Oz -o $output1.Ah_$seq.combined.vcf.gz
bcftools concat --allow-overlaps --threads $thr $output1.Dh_$seq.variant.recode.vcf.gz $output1.Dh_$seq.invariant.recode.vcf.gz -Oz -o $output1.Dh_$seq.combined.vcf.gz

parallel tabix {} ::: $output1.*h_$seq.*combined.vcf.gz