### ROH running 

#### prepare bi-allelic file for ROH caculation 
```
# use the output from Pixy filtered biallelic vcfs from each chromosome
# /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined_bi
module load picard/2.27.4
output=MKGDPRYuan_n121

picard GatherVcfs \
$(for vcf in /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined_bi/*.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O $output.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s ../rename.MKGDPRYuan_n121.txt $output.AhDh.combined.bi.vcf.gz -o $output.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" $output.AhDh.combined.bi.rehead.vcf >  $output.AhDh.combined.bi.rehead.id.vcf

grep -v '#' $output.AhDh.combined.bi.rehead.id.vcf | wc -l
vcftools --vcf $output.AhDh.combined.bi.rehead.id.vcf --het --out $output.het

```

#### setup vcf format for ROH running
```
Dir=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined_bi
vcf=MKGDPRYuan_n121.AhDh.combined.bi.rehead.id.vcf
output=MKGDPRYuan_n121

# rename the chromosome from At and Dt to 01 to 26
ml vcftools bcftools
grep ">" /lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa | sed 's/>//g' | sort | uniq | awk '{printf "%s\t%d\n", $0, NR}' > rename.chr.txt
bcftools annotate --rename-chrs rename.chr.txt $Dir/$vcf -o $output.AhDh.combined.bi.rehead.id.rech.vcf

# convert vcf into ped format via plink
module purge
module load plink/1.90b6.21
plink --vcf $output.AhDh.combined.bi.rehead.id.rech.vcf --allow-extra-chr --make-bed --const-fid --recode --out $output

module load r
Rscript ROH.R
```

#### Rcode below used the Rpackage [detectRUNS](https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html)
```
setwd(getwd())
library(detectRUNS)

slidingRuns <- slidingRUNS.run(
  genotypeFile = "MKGDPRYuan_n121.ped", 
  mapFile = "MKGDPRYuan_n121.map", 
  windowSize = 15,   #the size of sliding window (number of SNP loci) (default = 15)
  threshold = 0.05,  #the threshold of overlapping windows of the same state (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)
  minSNP = 10,       #minimum n. of SNP in a RUN (default = 3)
  ROHet = FALSE,     
  maxOppWindow = 1,  #max n. of homozygous/heterozygous SNP in the sliding window (default = 1)
  maxMissWindow = 1, #max. n. of missing SNP in the sliding window (default = 1)
  maxGap = 10^6,     #max distance between consecutive SNP to be still considered a potential run (default = 10^6 bps) 
  minLengthBps = 250000,   #minimum length of run in bps (defaults to 1000 bps = 1 kbps)
  minDensity = 1/10^3, #minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
  maxOppRun = NULL,
  maxMissRun = NULL
) 

# Runs of heterozygosity -- results were not plotted
slidingRuns_het <- slidingRUNS.run(
  genotypeFile = "MKGDPRYuan_n121.ped", 
  mapFile = "MKGDPRYuan_n121.map",
  windowSize = 10, 
  threshold = 0.05,
  minSNP = 10, 
  ROHet = TRUE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 10000, 
  minDensity = 1/10^6, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

save.image(file = "ROH_n86.final.RData")
```

#### For plotting details using RData see [Final_LD_plot.R](https://github.com/Wendellab/CaribbeanAD1/blob/main/07_ROH_He_121/Final_LD_plot.R)


#######################################################################

### Heterozygoisty caculating via vcftools
```
vcftools --vcf MKGDPRYuan_n121.AhDh.combined.bi.rehead.id.vcf --het --out MKGDPRYuan_n121.het
```
