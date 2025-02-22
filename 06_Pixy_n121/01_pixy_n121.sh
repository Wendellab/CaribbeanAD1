#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G 
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.pixy_n121.%J.out"
#SBATCH --job-name="pixy_n121"
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
seqA="Ah_"$seq
seqD="Dh_"$seq
thr=20 #NUMBER_THREADS
output1=MKGDPRYuan_n121
vcf=/work/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/00_VCFsubset/MKGDPRYuan_n121_combined/MKGDPRYuan_n121.AhDh.combined.rehead.vcf.gz 


######### random downsizing n = 5 for 20 replicates 
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


for i in {1..20}; do
    echo $i
    for sample in AD1_{MK,GD,PR_CR,PR_Ph,PR_PR325,Cultivar,LR1,LR2,Wild}; do
        grep "$sample" pixy_populationlist.txt | shuf -n 5 
    done > pixy_populationlist_$i.txt
done


######### caculate the pixy 


module purge
module load micromamba/1.4.2-7jjmfkf
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

mv $TMPDIR/$output1.*${seq}.subset${i}* /work/LAS/jfw-lab/weixuan/07_PRGD_popgene/04_MKPRGD/01_pixy/pixyoutput_n121/
done

micromamba deactivate
