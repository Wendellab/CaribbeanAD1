### Identifying the feral cottons in Caribbean Gh collections

#### Merging VCFs from selected samples

```
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
```

#### Final VCF for 145 samples 
```
module load picard/2.27.4

picard GatherVcfs \
$(for vcf in *.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O AD1AD2AD4_n145.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename.AD1AD2AD4_n145.txt AD1AD2AD4_n145.AhDh.combined.bi.vcf.gz -o AD1AD2AD4_n145.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" AD1AD2AD4_n145.AhDh.combined.bi.rehead.vcf >  AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf
```

#### Run Plink 
```
ml bcftools

module load plink/1.90b6.21

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/03_AD1AD2AD4_n145/AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf
output=AD1AD2AD4_n145

#Filtering LD by Plink 
plink --threads 30 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 30 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

plink --threads 30 --vcf $vcf --allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```

### A nice code block for R map plot of (Caribbean islands cotton samples)[https://github.com/Wendellab/CaribbeanAD1/blob/main/02_FeralCottonPR_n145/FINAL_GoogleMap_DistributionMapPlot.R]
