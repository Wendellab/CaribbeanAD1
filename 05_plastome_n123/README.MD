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
shopt -s globstar  ##globstar will let us match files recursively
files=( */*.fasta )  ##Array containing matched files, mention where to search and what files here

for i in "${files[@]}"; do d="${i%/*}"; d_="${d##*/}";  f="${i##*/}"; \
cp "$i" /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/plastome_output_aln"${d_}"_"$f"; done
```

#### renaming the plastome output fasta headers according to the file names and add "_1" and "_2" for two extraction results
```
for FILE in *.fasta;
do
 awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $FILE > changed_${FILE}
done

for FILE in changed*.fasta;
	do
	awk 'BEGIN{RS=">"}{if(NR>1)print ">"$1"_"(NR-1)"\n"$2}' $FILE > sorted_${FILE}
done
```

# some steps using mafft and eye scanning the two copies of cpDNA outcome not mentioned in here

###################################################################################

### For whole plastome

#### set up a dummpy IRB sequence from the annotation of AD1_GD_G7C via [Geseq](https://chlorobox.mpimp-golm.mpg.de/geseq.html), using mafft to and seqkit to remove IRB region across all samples, then remove the 'dummpy IRB' region sequence using sed 
```
module load mafft
mafft --thread 80 plastome12_IRB_n158.fasta > plastome12_IRB_n158_aln.fasta
seqkit mutate -d 89419:115020 --quiet -w 0  plastome12_IRB_n158_aln.fasta  | sed '/IRB/{N;d;}'  > plastome12_IRB_n158_aln_trimmed.fasta
```


#### building a nex format for MP based (PAUP4*)[https://paup.phylosolutions.com/]
```
#realign the sequences that without IRB region
mafft --thread 80 plastome12_IRB_n158_aln_trimmed.fasta > plastome12_IRB_n158_aln_trimmed_aln.fasta

#trim all gaps using trimal
/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/trimal/source/trimal -nogaps \
-in plastome12_IRB_n158_aln_trimmed_aln.fasta  -out plastome12_IRB_n158_aln_trimmed_aln_nogap.fasta

#sort the fasta based on sequence headers (names)
seqkit sort -n plastome12_IRB_n158_aln_trimmed_aln_nogap.fasta > plastome12_IRB_n158_aln_trimmed_aln_nogap_sort.fasta

#check the sequence duplication details to report in Haplotype analysis using [popart](https://popart.maths.otago.ac.nz/) 
cat plastome12_IRB_n158_aln_trimmed_aln_nogap_sort.fasta | seqkit rmdup -s -w 0 -o plastome_rmconst_rmdup.fasta -d plastome_rmconst_dup.fasta -D duplicated.detail.txt

#convert fasta into nex for PAUP and Popart
seqmagick convert --output-format nexus --alphabet dna plastome12_IRB_n158_aln_trimmed_aln_nogap_sort.fasta plastome12_IRB_n158_aln_trimmed_aln_nogap_sort.nex
```

#### ML based [IQtree2](https://github.com/iqtree/iqtree2)
```
module load iqtree2/2.2.2.7-py310-kcqit6n
iqtree2 -s ../plastome12_IRB_n158_aln_trimmed_aln_nogap_sort.fasta -B 1000 -T AUTO --prefix plastome12_noIRB_n158
```


###################################################################################

### For gene sequences in plastome

#### Transferring the annotation using AD1_GD_G7C results to all other cotton samples via [PGA](https://github.com/quxiaojian/PGA) 
```
# convert the alignment to unaligned fasta
sed 's/-//g' plastome12_n457_aln_sort.fasta > plastome12_n457_aln_sort_fixed.fasta
# seperate a multi-sample fasta into individual sequences
seqkit split --by-id --line-width 0 plastome12_n457_aln_sort_fixed.fasta

module load ncbi-rmblastn/2.14.0-py310-vqnew3z
module load blast

#Transfer the annotation from reference folder into targeted folder plastome457 and save the gb format output in putput PGA
perl /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/PGA/PGA.pl \
-r /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/PGA/plastomes/reference \
-t /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/PGA/plastomes/plastome457 \
-o /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/PGA/plastomes/outputPGA
```


#### extract the gene sequences using a customized python scrip [extract_gene_final.py](https://github.com/Wendellab/CaribbeanAD1/blob/main/05_plastome_n123/extract_gene_final.py) that takes a list of genebank format files and output the loci sequences in the output folder ./loci110_n158/
```
# generate a list of targeted gb files and extract the sequences
grep -f select_n158.txt plastome_n457.list > select_n158_gb.txt
python extract_gene_final.py select_n158_gb.txt ./loci110_n158/

```

#### build alignment for each gene sequence had >= 158 samples using mafft and remove the gaps
```
grep ">" -c *.fasta | sed 's/.fasta:/ /' | sort -g -k 2 | awk '$2 >= 158  { print $1}' > plastome_gene110.txt
word=$(cat plastome_gene110.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load mafft/7.508
echo "${word}"

mafft --auto --thread 10 "${word}".fasta > aln/"${word}"_aln.fasta
/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/trimal/source/trimal -nogaps -in "${word}"_aln.fasta -out "${word}"_trim.fasta  
```

#### Build the concatenated ML iqtree using the final trimmed alignments
```
module load iqtree2/2.2.2.7-py310-kcqit6n
iqtree2 -p ../aln --prefix cotton_n457_concat -B 1000 -T AUTO
```

#### concatenate all gene sequences into one fasta for PAUP and PopArt to run
```
seqkit concat ../aln/*.fasta > ./loci110_concat.fasta
```
