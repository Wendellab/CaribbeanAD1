### Selecting the Gb wild samples from [Yuan et al](https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202003634) samples 

#### Merging VCFs from selected samples

```
ml vcftools bcftools

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD5/AD5.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Ah_$seq.combined.vcf.gz \
-Oz -o MKAD2AD5AD4_n51.Ah_$seq.combined.vcf.gz

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD5/AD5.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Dh_$seq.combined.vcf.gz \
-Oz -o MKAD2AD5AD4_n51.Dh_$seq.combined.vcf.gz

bcftools view MKAD2AD5AD4_n51.Ah_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o MKAD2AD5AD4_n51.Ah_$seq.combined.bi.vcf.gz

bcftools view MKAD2AD5AD4_n51.Dh_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o MKAD2AD5AD4_n51.Dh_$seq.combined.bi.vcf.gz

module load parallel/20220522-sxcww47
parallel tabix {} ::: MKAD2AD5AD4_n51.*h_$seq.combined.bi.vcf.gz
```

#### Final VCF

```
module load picard/2.27.4

picard GatherVcfs \
$(for vcf in *.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename.MKAD2AD4AD5_n51.txt MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz -o MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf >  MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf

grep -v "#" MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf | wc -l
```


#### Run Plink to get PCA and ped file 

```
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
```


#### Run LEA for structure

```
module load r

Rscript LEA.R

cut -d ' ' -f 2 *.ped  > samplename.txt
```

