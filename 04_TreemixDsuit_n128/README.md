### Wild cotton relationships

#### Subsetting VCFs from selected samples
```
ml vcftools bcftools
output=AD1AD2AD4_n128

bcftools view -S subset_n128.txt /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/03_AD1AD2AD4_n145/AD1AD2AD4_n145.AhDh.combined.bi.rehead.id.vcf \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor \
-o $output.AhDh.combined.vcf
```


#### Purning LD using Plink
```
module load plink/1.90b6.21
vcf=$output.AhDh.combined.vcf

#Filtering LD by Plink 
plink --threads 10 --vcf $vcf  \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 10 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf
```


#### Setup for Treemix run
```
#mkdir treemix_GD2
cd treemix_GD2
mv /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/01_Plink_n128/$output.vcf ./

# convert vcf to Plink format
vcftools --vcf  $output.vcf --plink-tped --out  $output

# add the genetic group information in Plink tfam file, with all four GD2, and PR_Ph6 as seperate groups

awk '{gsub("AD1_GD_G1A", "AD1_GD2_G1A", $1); \
gsub("AD1_GD_G23A", "AD1_GD2_G23A", $1); \
gsub("AD1_GD_G24A", "AD1_GD2_G24A", $1); \
gsub("AD1_GD_G25B", "AD1_GD2_G25B", $1); \
gsub("AD1_PR_Ph_6", "AD1_PR_Ph2_6", $1); \
print}' $output.tfam | \
awk '{
    if ($1 ~ /^AD1_PR/) {
        match($1, /^[^_]*_[^_]*_[^_]*/)
        printf substr($1, RSTART, RLENGTH)
    } else {
        match($1, /^[^_]*_[^_]*/)
        printf substr($1, RSTART, RLENGTH)
    }
    for (i = 2; i <= NF; i++) {
        printf "\t%s", $i
    }
    printf "\n"
}' >  $output.2.tfam

# Use the new tfam file to generate a population list for treemix
cat $output.2.tfam |awk '{print $1"\t"$2"\t"$1}' > sample.pop.cov

# also change the file format for tped file to match iwth tfam file
cp $output.tped $output.2.tped

# caculate the allele frequency for each population/group using plink and gzip the outcome
plink --threads 10 --tfile $output.2 --freq --allow-extra-chr --within sample.pop.cov 
gzip plink.frq.strat

# convert plink output format into treemix format using treemix python code
python /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/plink2treemix.py plink.frq.strat.gz sample.treemix.in.gz
```




#### Running Treemix
```
mkdir Treemixoutput

module purge
module load gsl/2.7.1-uuykddp
module load python/3.10.10-zwlkg4l
module load boost/1.81.0-zwxu2hi
module load py-numpy/1.26.3-py310-gntgk2n

# m = migration events,
# and i = number of replicates per m -k
# *-se*  tells TreeMix to use standard errors when estimating covariance matrix of allele frequencies (recommended for bootstrapping).
# Sets *AD4_mus* as the root of the tree. This should be an outgroup population.
# *-k* Block size in SNPs for accounting for linkage disequilibrium. 1000 is typical, but depends on genome structure and SNP density.
for m in {0..8}
do
    {
    for i in {1..5}
    do
        {
            /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/src/treemix -se -bootstrap \
			-i sample.treemix.in.gz -root AD4_mus \
			-o Treemixoutput/TreeMix.${m}.${i} \
			-m ${m} -k 1000 
        }
    done
    }
done
```

#### select the best treemix model using OptM
```
#below are the codes in R
library(OptM)
linear = optM("./Treemixoutput")
plot_optM(linear, plot = F, pdf = "file.pdf")
```

#### Setup for Dsuite
```
mkdir Dsuite

# using Treemix to get a ML tree for Dsuite
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/src/treemix -i sample.treemix.in.gz -o sample.ML.tree -root AD4_mus -k 1000
zcat sample.ML.tree.treeout.gz | sed 's/AD4_mus/Outgroup/g' > Dsuite/sp.tre

# using popinformation from treemix to generate a individual list for Dsuite
cut -f 2,3 sample.pop.cov | awk '{gsub("AD4_mus", "Outgroup", $2); print}'  > Dsuite/sets.txt
sed -i 's/ /\t/g' Dsuite/sets.txt
```

#### Running Dsuite
```
cd Dsuite

/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/Build/Dsuite Dtrios -c -n $output -t sp.tre ../$output.vcf sets.txt 
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/Build/Dsuite  Fbranch  sp.tre sets_${output}_tree.txt > species_sets_${output}_Fbranch.txt
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/utils/dtools.py species_sets_${output}_Fbranch.txt sp.tre --dpi 300 --tree-label-size 18 --ladderize
```
