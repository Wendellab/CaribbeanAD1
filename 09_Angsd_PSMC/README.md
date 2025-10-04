### Tajima's D and Watterson's theta estimation using [Angsd Thetas](https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests) and [Angsd SFS](https://www.popgen.dk/angsd/index.php/SFS_Estimation)

#### MAF SFS and Tajima's D and Watterson's theta
```
module load angsd
cd $DIR
angsd -bam bamlist_$file.txt -out $name -anc $ref -P $thr -doSaf 1 -doMaf 1 -doMajorMinor 1 -doGlf 3 -uniqueOnly -GL 2 -minMapQ 30 -minQ 20 -minInd 15

#-doMaf 1: Frequency (fixed major and minor)
#-doMajorMinor 1: Infer major and minor from GL
#-uniqueOnly: reads mapped to unique location
#GL 2: genotype likelihood model
#mapped quality filtering -minMapQ 30 -minQ 20 -minInd 25
#-doSaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
#-doGlf 3: binary beagle in/out-put gz requires -doMajorMinor [1|2]


realSFS $name.saf.idx -maxIter 100 -P $thr -fold 1 > $name.sfs
realSFS saf2theta $name.saf.idx -outname $name -sfs $name.sfs -fold 1
thetaStat do_stat $name.thetas.idx
thetaStat do_stat $name.thetas.idx -win 50000 -step 10000  -outnames $name.thetasWindow.gz
```

### Ne estimation using [PSMC](https://github.com/lh3/psmc)
#### Setup the format from bam to PSMC input format
```
ref=/work/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=30 #NUMBER_THREADS
outdir=01_psmcoutput

file=$(cat bamfile_list.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file .realign.bam)

module load vcftools
module load bcftools
module load psmc
module load samtools/1.16.1

module load psmc/2016-1-21-ln6nqee
module load gnuplot/5.4.3-py310-wirgelt
module load ghostscript

samtools mpileup -C50 -uf $ref $file | bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > $name.fq.gz
#Creates pileup from BAM file ($file) using reference genome ($ref)
#-C50: Adjusts mapping quality for reads with excessive mismatches (for better indel handling)
#-u: Uncompressed BCF output
#-f: Reference genome

#Calls variants in consensus mode (-c)

#Converts VCF to a consensus fastq file with: Minimum depth 10 (-d 10), Maximum depth 100 (-D 100) — filters out likely PCR duplicates or repetitive regions

fq2psmcfa -q20 $name.fq.gz > $name.psmcfa
#Converts the fastq to a psmcfa file, which is a pseudo-diploid fasta-like format required by PSMC
#-q20: Filters out bases with quality scores below 20

splitfa $name.psmcfa > split.$name.psmcfa
#Splits the input into non-overlapping segments (typically 100) for bootstrapping
#Bootstrapping is used to assess confidence intervals in PSMC estimates
```

#### Running PSMC with replicates [PSMC](https://github.com/lh3/psmc)
```
#SBATCH --array=1-100

ref=/work/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=30 #NUMBER_THREADS
outdir=01_psmcoutput
seq=${SLURM_ARRAY_TASK_ID}

module load vcftools
module load bcftools
module load psmc
module load samtools/1.16.1

module load psmc/2016-1-21-ln6nqee
module load gnuplot/5.4.3-py310-wirgelt
module load ghostscript

psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" split.CR_1.psmcfa -o $outdir/round${seq}.CR_1.psmc
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" split.Ph_3.psmcfa -o $outdir/round${seq}.Ph_3.psmc
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" split.Pop1_1.psmcfa -o $outdir/round${seq}.Pop1_1.psmc
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" split.PR325_15.psmcfa -o $outdir/round${seq}.PR325_15.psmc
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" split.YUC-G17A.psmcfa -o $outdir/round${seq}.YUC-G17A.psmc
```

### Ne estimation using [SMC++](https://github.com/popgenmethods/smcpp)
#### Prepar input VCF file for smc++
```
ml vcftools bcftools

output=MKGDPR_n81

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined_bi/MKGDPRYuan_n121.AhDh.combined.bi.rehead.vcf

bcftools query -l $vcf | grep -E "MK|GD|PR" > subset_n81.txt

bcftools view -S subset_n81.txt $vcf \
-m2 -M2 -i 'F_MISSING==0' -q 0.001:minor \
-o $output.AhDh.combined.bi.rehead.vcf

bgzip $output.AhDh.combined.bi.rehead.vcf
tabix $output.AhDh.combined.bi.rehead.vcf.gz

bgzip AD1.TX2094.v2.masked.bed
tabix -p bed AD1.TX2094.v2.masked.bed.gz
```


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G 
#SBATCH --time=05:30:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.beagle_flare.%J.out"
#SBATCH --job-name="beagle_flare"
#SBATCH --array=1-13

ml vcftools bcftools

chr=$(printf %02d ${SLURM_ARRAY_TASK_ID})
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/06_smcpp/MKGDPR_n81.AhDh.combined.bi.rehead.vcf.gz

echo Ah_$chr
echo Dh_$chr

module purge
module load openjdk/21.0.3_9-vngib7s

java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$vcf out=MKGDPR_n81.Ah_$chr.phased chrom=Ah_$chr
java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$vcf out=MKGDPR_n81.Dh_$chr.phased chrom=Dh_$chr

tabix MKGDPR_n81.Ah_$chr.phased.vcf.gz
tabix MKGDPR_n81.Dh_$chr.phased.vcf.gz

bcftools concat -Oz -o ../MKGDPR_n81.AhDh.combined.phased.vcf.gz $(for vcf in *phased*.gz; do echo "$vcf"; done)

cd ../

# 1. Extract contig lengths from the reference VCF
bcftools view -h MKGDPR_n81.AhDh.combined.bi.rehead.vcf.gz | \
grep '^##contig' | \
awk -F'[=<,]' '{for(i=1;i<=NF;i++){if($i=="ID"){id=$(i+1)}; if($i=="length"){len=$(i+1)}}; if(id && len){print id "\t" len; id=""; len=""}}' > contig_lengths.txt

# contig_lengths.txt will look like:
# Ah_01   119153486
# Ah_02   108577906
# ...

# 2. Replace contig lines in the target VCF
while IFS=$'\t' read -r contig length; do
    sed -i "s/##contig=<ID=$contig[^>]*>/##contig=<ID=$contig,length=$length,assembly=unknown>/" MKGDPR_n81.AhDh.combined.phased.vcf
done < contig_lengths.txt

bgzip MKGDPR_n81.AhDh.combined.phased.vcf
tabix -f MKGDPR_n81.AhDh.combined.phased.vcf.gz
```

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G 
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.smcpp.%J.out"
#SBATCH --job-name="smcpp"

####
module load bcftools
vcf=MKGDPR_n81.AhDh.combined.phased.vcf.gz

MKlist=$(bcftools query -l $vcf | grep 'MK' | paste -sd,)
GDlist=$(bcftools query -l $vcf | grep 'GD' | paste -sd,)
PR_CRlist=$(bcftools query -l $vcf | grep 'PR_CR' | paste -sd,)
PR_Phlist=$(bcftools query -l $vcf | grep 'PR_Ph' | paste -sd,)
PR_PR325list=$(bcftools query -l $vcf | grep 'PR_PR325' | paste -sd,)

echo $MKlist
echo $GDlist
echo $PR_CRlist
echo $PR_PR325list
echo $PR_Phlist
module purge

####
source /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/06_smcpp/smcpp_env/bin/activate 

module load gsl/2.7.1-uuykddp
module load mpfr/4.2.0-av6zkh3
module load gmp/6.2.1-dh2y7b4
module load gcc/12.2.0-khmr45w
module load python/3.10.10-zwlkg4l
module load py-numpy/1.26.3-py310-gntgk2n

####
mkdir -p 01_output

for i in {01..13}
do
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MK.Ah_${i}.smc.gz Ah_${i} MK:$MKlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MK.Dh_${i}.smc.gz Dh_${i} MK:$MKlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GD.Ah_${i}.smc.gz Ah_${i} GD:$GDlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GD.Dh_${i}.smc.gz Dh_${i} GD:$GDlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/CR.Ah_${i}.smc.gz Ah_${i} CR:$PR_CRlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/CR.Dh_${i}.smc.gz Dh_${i} CR:$PR_CRlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/Ph.Ah_${i}.smc.gz Ah_${i} Ph:$PR_Phlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/Ph.Dh_${i}.smc.gz Dh_${i} Ph:$PR_Phlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/PR325.Ah_${i}.smc.gz Ah_${i} PR325:$PR_PR325list
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/PR325.Dh_${i}.smc.gz Dh_${i} PR325:$PR_PR325list
done


# Set mutation rate
MUT_RATE=4.56e-9

# Create analysis directories
mkdir -p 02_analysis/MK
mkdir -p 02_analysis/GD
mkdir -p 02_analysis/CR
mkdir -p 02_analysis/Ph
mkdir -p 02_analysis/PR325

# Run smc++ estimate for each population in parallel
smc++ estimate --timepoints 10 1000000 --cores 30 -o 02_analysis/MK $MUT_RATE 01_output/MK.*.smc.gz 
smc++ estimate --timepoints 10 1000000 --cores 30 -o 02_analysis/GD $MUT_RATE 01_output/GD.*.smc.gz 
smc++ estimate --timepoints 10 1000000 --cores 30 -o 02_analysis/CR $MUT_RATE 01_output/CR.*.smc.gz 
smc++ estimate --timepoints 10 1000000 --cores 30 -o 02_analysis/Ph $MUT_RATE 01_output/Ph.*.smc.gz 
smc++ estimate --timepoints 10 1000000 --cores 30 -o 02_analysis/PR325 $MUT_RATE 01_output/PR325.*.smc.gz 

###
cd 02_analysis/
smc++ plot plot.pdf -g 1 */*.final.json --ylim 0 1000000 -c

###

```


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G 
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.smcpp.%J.out"
#SBATCH --job-name="smcpp"

####
module load bcftools
vcf=MKGDPR_n81.AhDh.combined.phased.vcf.gz

MKlist=$(bcftools query -l $vcf | grep 'MK' | paste -sd,)
GDlist=$(bcftools query -l $vcf | grep 'GD' | paste -sd,)
PR_CRlist=$(bcftools query -l $vcf | grep 'PR_CR' | paste -sd,)
PR_Phlist=$(bcftools query -l $vcf | grep 'PR_Ph' | paste -sd,)
PR_PR325list=$(bcftools query -l $vcf | grep 'PR_PR325' | paste -sd,)

echo $MKlist
echo $GDlist
echo $PR_CRlist
echo $PR_PR325list
echo $PR_Phlist
module purge

####
source /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/06_smcpp/smcpp_env/bin/activate 

module load gsl/2.7.1-uuykddp
module load mpfr/4.2.0-av6zkh3
module load gmp/6.2.1-dh2y7b4
module load gcc/12.2.0-khmr45w
module load python/3.10.10-zwlkg4l
module load py-numpy/1.26.3-py310-gntgk2n

####
mkdir -p 01_output/MKGD
mkdir -p 01_output/MKPh
mkdir -p 01_output/GDPh

for i in {01..13}
do
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKGD/MK_GD.Ah_${i}.smc.gz Ah_${i} MK:$MKlist GD:$GDlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKGD/MK_GD.Dh_${i}.smc.gz Dh_${i} MK:$MKlist GD:$GDlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKGD/GD_MK.Ah_${i}.smc.gz Ah_${i} GD:$GDlist MK:$MKlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKGD/GD_MK.Dh_${i}.smc.gz Dh_${i} GD:$GDlist MK:$MKlist


smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKPh/Ph_MK.Ah_${i}.smc.gz Ah_${i} Ph:$PR_Phlist MK:$MKlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKPh/Ph_MK.Dh_${i}.smc.gz Dh_${i} Ph:$PR_Phlist MK:$MKlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKPh/MK_Ph.Ah_${i}.smc.gz Ah_${i} MK:$MKlist Ph:$PR_Phlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/MKPh/MK_Ph.Dh_${i}.smc.gz Dh_${i} MK:$MKlist Ph:$PR_Phlist


smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GDPh/Ph_GD.Ah_${i}.smc.gz Ah_${i} Ph:$PR_Phlist GD:$GDlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GDPh/Ph_GD.Dh_${i}.smc.gz Dh_${i} Ph:$PR_Phlist GD:$GDlist

smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GDPh/GD_Ph.Ah_${i}.smc.gz Ah_${i} GD:$GDlist Ph:$PR_Phlist
smc++ vcf2smc --mask AD1.TX2094.v2.masked.bed.gz $vcf 01_output/GDPh/GD_Ph.Dh_${i}.smc.gz Dh_${i} GD:$GDlist Ph:$PR_Phlist

done

mkdir -p 02_analysis/MKGD
mkdir -p 02_analysis/MKPh
mkdir -p 02_analysis/GDPh

smc++ split --timepoints 10 1000000 --cores 30 -o 02_analysis/MKGD/ 02_analysis/MK/model.final.json 02_analysis/GD/model.final.json 01_output/MKGD/*.smc.gz
smc++ split --timepoints 10 1000000 --cores 30 -o 02_analysis/MKPh/ 02_analysis/MK/model.final.json 02_analysis/Ph/model.final.json 01_output/MKPh/*.smc.gz
smc++ split --timepoints 10 1000000 --cores 30 -o 02_analysis/GDPh/ 02_analysis/GD/model.final.json 02_analysis/Ph/model.final.json 01_output/GDPh/*.smc.gz

 
```
