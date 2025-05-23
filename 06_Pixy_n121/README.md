### Pixy running with downsampling for Caribbean cottons

#### Generate 20 lists with each population downsample for 5 individuals 
```
# extract the first two parts of strings before second underscore as population name
bcftools query -l \
/work/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined/MKGDPRYuan_n121.AhDh.combined.rehead.vcf.gz | \
awk '{
    if ($1 ~ /^AD1_PR/) {
        match($1, /^[^_]*_[^_]*_[^_]*/)
        prefix = substr($1, RSTART, RLENGTH)
    } else {
        match($1, /^[^_]*_[^_]*/)
        prefix = substr($1, RSTART, RLENGTH)
    }
    printf "%s\t%s", $1, prefix
    for (i = 2; i <= NF; i++) {
        printf "\t%s", $i
    }
    printf "\n"
}' > pixy_populationlist.txt

# downsize for 5 individual for 20 replicates
for i in {1..20}; do
    echo $i
    for sample in AD1_{MK,GD,PR_CR,PR_Ph,PR_PR325,Cultivar,LR1,LR2,Wild}; do
        grep "$sample" pixy_populationlist.txt | shuf -n 5 
    done > pixy_populationlist_$i.txt
done
```

#### Run Pixy chromsome 
```
seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
seqA="Ah_"$seq
seqD="Dh_"$seq
thr=20 #NUMBER_THREADS
output1=MKGDPRYuan_n121
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined/MKGDPRYuan_n121.AhDh.combined.rehead.vcf.gz 

module purge
module load micromamba/1.4.2-7jjmfkf
#export MAMBA_ROOT_PREFIX=/home/weixuan
eval "$(micromamba shell hook --shell=bash)"
micromamba activate

for i in {1..20}
do
echo $i
echo "Ah_"$seq
echo "Dh_"$seq
echo $vcf
pixy --stats pi fst dxy --bypass_invariant_check yes --chromosomes $seqA --vcf $vcf --populations pixy_populationlist_$i.txt --window_size 10000 --n_cores $thr --output_folder $TMPDIR/ --output_prefix $output1.Ah$seq.subset$i
pixy --stats pi fst dxy --bypass_invariant_check yes --chromosomes $seqD --vcf $vcf --populations pixy_populationlist_$i.txt --window_size 10000 --n_cores $thr --output_folder $TMPDIR/ --output_prefix $output1.Dh$seq.subset$i

mv $TMPDIR/$output1.*${seq}.subset${i}* /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/01_pixy/pixyoutput_n121/
done

micromamba deactivate
```

#### Caculate the average Pixy for each window accross all replicates for dxy, pi and fst, respectively 
```
awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_dxy.txt > output_temp/MKGDPRYuan_n121.dxy.txt
awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_pi.txt > output_temp/MKGDPRYuan_n121.pi.txt
awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_fst.txt > output_temp/MKGDPRYuan_n121.fst.txt

cd output_temp

cut -f 1,2,5 MKGDPRYuan_n121.pi.txt > MKGDPRYuan_n121.pi2.txt
cut -f 1,2,3,6 MKGDPRYuan_n121.dxy.txt > MKGDPRYuan_n121.dxy2.txt
cut -f 1,2,3,6 MKGDPRYuan_n121.fst.txt > MKGDPRYuan_n121.fst2.txt

awk 'FNR > 1 {
    group = $1"\t"$2        # Create a key for the group using the first and second columns
    count[group]++          # Increment the count for this group
    sum[group] += $3        # Accumulate the sum of the third column for this group
    sumsq[group] += $3*$3   # Accumulate the square of the third column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                     # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0             # Calculate standard deviation
        print group, mean, sd                                # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.pi2.txt | sort  > MKGDPRYuan_n121.pi3.txt


awk 'FNR > 1 {
    group = $1"\t"$2"\t"$3   # Create a key for the group using the first, second, and third columns
    count[group]++           # Increment the count for this group
    sum[group] += $4         # Accumulate the sum of the fourth column for this group
    sumsq[group] += $4*$4    # Accumulate the square of the fourth column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                      # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0              # Calculate standard deviation
        print group, mean, sd                                 # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.dxy2.txt | sort  > MKGDPRYuan_n121.dxy3.txt

awk 'FNR > 1 {
    group = $1"\t"$2"\t"$3   # Create a key for the group using the first, second, and third columns
    count[group]++           # Increment the count for this group
    sum[group] += $4         # Accumulate the sum of the fourth column for this group
    sumsq[group] += $4*$4    # Accumulate the square of the fourth column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                      # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0              # Calculate standard deviation
        print group, mean, sd                                 # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.fst2.txt | sort  > MKGDPRYuan_n121.fst3.txt
```

