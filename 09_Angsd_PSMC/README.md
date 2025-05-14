### Tajima's D and Watterson's theta estimation using [Angsd Thetas] (https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests) and [Angsd SFS](https://www.popgen.dk/angsd/index.php/SFS_Estimation)

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
### Setup the format from bam to PSMC input format
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


