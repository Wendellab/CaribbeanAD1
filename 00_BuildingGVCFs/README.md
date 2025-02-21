### A353 loci sequence extraction and gene trees reconstruction
#### Raw reads trimming 
```
module load Trimmomatic/0.39-Java-1.8.0_144

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/01_trimmed/*_1.fq.gz
	do
	withpath="${file}"
	filename=${withpath##*/}
	base="${filename%**_1.fq.gz}"
	echo "${base}"
	trimmomatic PE -threads 25 "${base}"_1.fq.gz "${base}"_2.fq.gz \
	"${base}"_1P.fq.gz "${base}"_1UP.fq.gz \
	"${base}"_2P.fq.gz "${base}"_2UP.fq.gz \
	ILLUMINACLIP:TruSeq3-PE-2.fa:2:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done
```
