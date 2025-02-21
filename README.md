# Comparative population genomics of relictual Caribbean island *Gossypium hirsutum*

Manuscript link: to be updated

<p float="left">
  <img src="Supplementary/Cabo Rojo flower.JPEG" height="170" />
  <img src="Supplementary/Cabo Rojo population typical view.JPEG" height="170" /> 
  <img src="Supplementary/Cabo Rojo seeds and fibers 1.JPEG" height="170" /> 
  <img src="Supplementary/PR325 Salinas Providencia Population 11.JPEG" height="170" /> 
  <br>
	<img src="Supplementary/YUC_boll.jpg" height="170" />
	<img src="Supplementary/YUC_buds copy.jpg" height="170" /> 
	<img src="Supplementary/YUC_flower.jpg" height="170" /> 
</p>

#
### Abstract 
*Gossypium hirsutum* is the world’s most important source of cotton fiber, yet the diversity and population structure of its wild forms remain largely unexplored. The complex domestication history of G. hirsutum combined with its reciprocal introgression with a second domesticated species, *G. barbadense*, has generated a wealth of morphological forms and feral derivatives of both species and their interspecies recombinants, which collectively are scattered across a large geographic range in arid regions of the Caribbean basin. In this study, we aimed to assess genetic diversity within and among populations from two Caribbean islands, Puerto Rico (n = 43, five sites) and Guadeloupe (n = 25, one site), which contain putative wild and introgressed forms. Using whole genome resequencing data combined with a phylogenomic framework derived from a broader genomic survey, we parsed individuals into feral derivatives and truly wild forms. Feral cottons variously show genetic and morphological resemblance to the domesticated cottons, and vary greatly in genetic variation and heterozygosity, reflecting a complex history of interspecific and intraspecific gene flow that is spatially highly variable in its effects. Furthermore, wild cottons in both Caribbean islands appear to be relatively inbred, especially the Guadeloupe samples. Our results highlight the dynamics of population demographics in relictual wild cottons that experienced profound genetic bottlenecks associated with habitat destruction superimposed on a natural pattern of widely scattered populations. These results have implications for conservation of wild diversity in *G. hirsutum*. 

#
### Gene tree based  phylogeny inference 

We extracted Angiosperm353 loci via [HybPiper](https://github.com/mossmatters/HybPiper) using target-enrichment sequencing, and high-copy marker of plastome and nrDNA were extracted via [Getorgenalla](https://github.com/Kinggerm/GetOrganelle). All bioinformatic were performed via New Zealand eScience Infrastructure [NeSI](https://www.nesi.org.nz/)
 
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

#### Supercontigs retreving via HybPiper
```
module load Python/3.9.9-gimkl-2020a
module load HybPiper/2.0.1rc-Miniconda3
module load Parallel/20200522

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/*_1P.fq.gz
	do
	withpath="${file}"
	filename=${withpath##*/}
	base="${filename%*_1P.fq.gz}" 
	echo "${base}"
	
	hybpiper assemble --run_intronerate \
	--readfiles /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/"${base}"*P.fq.gz \
	--targetfile_dna mega353.fasta --bwa \
	--prefix "${base}"_nomerge  --no_padding_supercontigs \
	--timeout_assemble 600 --paralog_min_length_percentage 0.5
	done
```
### High copy gene sequence extraction
#### Raw reads trimming
```
module load Trimmomatic/0.39-Java-1.8.0_144

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/01_trimmedreads/*_1.fq.gz
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
#### Retreiving high-copy markers via Getorganelle 
```
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /scale_wlg_persistent/filesets/project/massey02696/biopython
cd /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/03_getorgenlle


for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_GenomeSkimming_NCBI_upload/02_trimmedreads/*_1P.fq.gz 
	do 
	withpath="${file}" 
	filename=${withpath##*/} 
	base="${filename%*_1P.fq.gz}" 
	echo "${base}" 
	
	get_organelle_from_reads.py -1 ../02_trimmedreads/"${base}"_1P.fq.gz -2 ../02_trimmedreads/"${base}"_2P.fq.gz -o plastome_output/"${base}" -R 15 -k 21,45,65,85,105 -F embplant_pt
	get_organelle_from_reads.py -1 ../02_trimmedreads/"${base}"_1P.fq.gz -2 ../02_trimmedreads/"${base}"_2P.fq.gz -o nrdna_output/"${base}" -R 10 -k 35,85,115 -F embplant_nr
	done
```
#
### SNPs based phylogeny analysis 
SNP amount Hyb-Seq reads were extracted using [HybSeq-SNP-Extraction](https://github.com/lindsawi/HybSeq-SNP-Extraction)

#### SNP calling via remapping the reads back to the same reference for 23 representative individuals
```
#The reference fasta gene names need to be changed according to file names
#The script variantcall need to change fastaq to fq
#The GenotypetoPCA need to change the script expression & to || for GATK
#The plink needs to update the version 

module rest
module load GATK/4.1.4.1-gimkl-2018b
module load PLINK/1.09b6.16
module load SAMtools/1.8-gimkl-2018b
module load BWA/0.7.17-gimkl-2017a
module load BCFtools/1.9-GCC-7.4.0
module load Python/3.7.3-gimkl-2018b
module load WhatsHap/1.1-gimkl-2020a

bash variantcall.sh Azho-AK16_S7_L001.supercontigs.fasta Azho-AK16_S7_L001 #please note that sample name do not have the slash / at the end
bash GenotypesToPCA.sh Azho-AK16_S7_L001.supercontigs.fasta Azho-AK16_S7_L001
bash plink_stats.sh Azho-AK16_S7_L001
```



