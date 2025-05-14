### Wild cotton plastome relationships

#### extracting whole plastomes using trimmed WGS reads
```
module load getorganelle/1.7.7.0-py310-u45ybv3

DIR=/lustre/hdd/LAS/jfw-lab/weixuan/03_Ggvcf/AD2_Wild/00_rawreads/trimmedReads
outputDIR=/lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle
thr=10

file=$(ls -1 $DIR/*.R1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file .R1.fq.gz)

# adjust -R and -k flags if needed. for example -R 200, and -k 125 worked for my samples 
get_organelle_from_reads.py -t $thr -1 $DIR/"${name}".R1.fq.gz -2 $DIR/"${name}".R2.fq.gz -o $outputDIR/plastome_output/"${name}" -R 15 -k 21,45,65,85,105 -F embplant_pt
```

#### renaming the plastome output file names by adding sample infor
```
for i in "${files[@]}"; do d="${i%/*}"; d_="${d##*/}";  f="${i##*/}"; \
cp "$i" /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/plastome_output_aln"${d_}"_"$f"; done
```

#### renaming the plastome output file names by adding sample infor
```
