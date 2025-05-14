### LDdecay analysis used [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay)

#### since the analysis is sensitive for minor allele frequency filtering and LD pruning, so we did not performed any additional filtering other than biallelic and depth, however, the sample size is small for this analysis. We also performed a downsampling (n=5) for this analaysis for 20 replicates
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined_bi/MKGDPRYuan_n121.AhDh.combined.bi.rehead.id.vcf

# get the subpopulation sample list
ml bcftools
bcftools query -l $vcf | grep 'MK' > MK.list
bcftools query -l $vcf | grep 'PR_Ph' > PR_Ph.list
bcftools query -l $vcf | grep 'PR_PR325' > PR_PR325.list
bcftools query -l $vcf | grep 'PR_CR' > PR_CR.list
bcftools query -l $vcf | grep 'Cultivar' > Cultivar.list
bcftools query -l $vcf | grep 'Wild' > Wild.list
bcftools query -l $vcf | grep 'LR1' > LR1.list
bcftools query -l $vcf | grep 'LR2' > LR2.list
bcftools query -l $vcf | grep 'GD' > GD.list
```

#### Running PR_Ph population with all n =5 
```
/lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat PR_Ph.stat.gz -SubPop PR_Ph.list -MaxDist 500
```


#### Downsize for all remaining populations
```
for i in {1..20}; do
  shuf -n 5 MK.list > MK.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat MK.list_"$i".stat.gz -SubPop MK.list_"$i" -MaxDist 500
  shuf -n 5 GD.list > GD.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat GD.list_"$i".stat.gz -SubPop GD.list_"$i" -MaxDist 500
  shuf -n 5 Wild.list > Wild.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat Wild.list_"$i".stat.gz -SubPop Wild.list_"$i" -MaxDist 500
  shuf -n 5 Cultivar.list > Cultivar.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat Cultivar.list_"$i".stat.gz -SubPop Cultivar.list_"$i" -MaxDist 500
  
  shuf -n 5 LR1.list > LR1.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat LR1.list_"$i".stat.gz -SubPop LR1.list_"$i" -MaxDist 500
  shuf -n 5 LR2.list > LR2.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat LR2.list_"$i".stat.gz -SubPop LR2.list_"$i" -MaxDist 500
  shuf -n 5 PR_CR.list > PR_CR.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat PR_CR.list_"$i".stat.gz -SubPop PR_CR.list_"$i" -MaxDist 500
  shuf -n 5 PR_PR325.list > PR_PR325.list_"$i"
  /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat PR_PR325.list_"$i".stat.gz -SubPop PR_PR325.list_"$i" -MaxDist 500

done
```

#### Tight up the format for plotting
```
# generate a list that contain files of LDdecay output
ls $PWD/*stat.gz | awk '{print $1, $1}' | awk '{gsub("/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/02_LD_decay/","", $2); gsub (".stat.gz", "", $2); print}' > Pop.ResultPath.list

# use PopLDdecay perl script to make the raw plot and save the output for R plot
module load r
perl /lustre/hdd/LAS/jfw-lab/weixuan/04_MK_Known_gvcf/VCF_output3/MK_Known_n65/LD_decay/PopLDdecay/bin/Plot_MultiPop.pl  -inList  Pop.ResultPath.list  -output output_n121/Fig -keepR


## this is write an additional line in the file by adding file names to merge all files together to plot in R
for file in Fig.*; do
    # Extract the base name of the file (without the extension)
    base_name=$(basename "$file")
    
    # Process the file and add the base name as a new column
    awk -v name="$base_name" 'NR==1 {print $0, "name_file"} NR>1 {print $0, name}' "$file" >> temp/LD_$base_name
done

awk 'FNR>1 || NR==1' LD_Fig.*list_* | sed 's/#//g; s/Fig.//g'  > LD_Fig_test.txt
```

#### For final plotting see R code [Final_LD_plot.R](https://github.com/Wendellab/CaribbeanAD1/blob/main/07_ROH_He_121/Final_LD_plot.R)

