---
output:
  html_document: default
  pdf_document: default
---
# Manual epiGBS 2.0

- [Prerequisites for running the pipeline](#prerequisites-for-running-the-pipeline)
- [Preparation to run the pipeline](#preparation-to-run-the-pipeline)
- [Start the pipeline](#start-the-pipeline)
- [Explanation of files in the output directory](#explanation-of-files-in-the-output-directory)
- [When not to run the pipeline?](#when-not-to-run-the-pipeline)
- [Quality control or "How to dicover errors?"](#quality-control-or-how-to-discover-errors)
- [Example Config Files](#example-config-files)
- [List of used software and references](#list-of-used-software-and-references)

## Prerequisites from the weblab for running the pipeline

- A basic knowledge of Linux
- Working with files and directories (cd, ls, nano)
- Executing commandos (git, conda, snakemake)
- Adapters include
	- Control nucleotide
	- Wobble (running without is also possible)
- No reference genome required
- Sequencing: paired-end Illumina-sequencing
- Readfiles: un-demultiplexing and un-trimmed
- Linux server (e.g. Ubuntu 16.04, Sudo rights not necessary)
- Miniconda
	- Download with `wget` https://docs.conda.io/en/latest/miniconda.html
	- Installing like described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), can be installed in user's home without sudo
		- `bash Miniconda3-latest-Linux-x86_64.sh` 

## Preparation to run the pipeline

- Make a conda environment for snakemake if snakemake is not installed globally on the server. You do not need administrator rights to do this but conda has to be installed. 
	- `conda create -n snake snakemake=5.4.5`
	- `conda activate snakemake`
- Make a copy of the pipeline
	- `git clone https://gitlab.bioinf.nioo.knaw.nl/FleurG/epigbs-snakemake.git epiGBS2.0`
- Enter the created directory:
	- `cd epiGBS2.0`
- Open and fill in the config file: __All paths are full paths, no relative paths allowed.__ 
	- `nano config.yaml`
	- Output_dir: Path of directory to store all output files and directory. Path will be created by the pipeline if it does not exist.
	- Input_dir: Path of directory containing raw data (e.g. fastq.gz, or fq.gz) and barcode file
	- Barcodes: Filename of the barcode file
	- Tmpdir: Path to the directory where temporary files will be stored. For most systems this will be /tmp
	- Cycles: read length
	- Read1: Name of the read file containing forward (R1) reads
	- Read2: Name of the read file containing forward (R2) reads
	- Ref_dir: Only needed for reference-based analysis. Path to directory containing a reference genome
	- Genome: Only needed for reference-based analysis. Name of the genome file (prefix .fa)
	- Mode: choice between ???denovo??? or ???reference???
	- Param_denovo: Optional and only relevant for analysis in de-novo mode. To run on default choose ?????? or ???default???, otherwise enter a number
	- Identity: 3rd clustering step, default: ???0.95???
	- Min-depth: First cluster-step. Minimal number of reads to build a cluster.
	- Max-depth: First cluster-step. Maximal number of reads to build a cluster.

- Make a barcode file: The barcode file is tab-delimited and contains at least the following columns: Flowcell, Lane, Barcode_R1, Barcode_R2, Sample, ENZ_R1, ENZ_R2, Wobble_R1, Wobble_R2. All other fields are optional. If you make the barcode file in Excel, make sure that no `^M` are present after uploading the file to the Linux server. You can check this by opening the barcode file with `cat -A barcode.tsv`.
```
# barcodes.tsv
Flowcell        Lane    Barcode_R1      Barcode_R2      Sample  history Country PlateName       Row     Column  ENZ_R1  ENZ_R2  Wobble_R1       Wobble_R2       Species
H53KHCCXY       5       AACT    CCAG    BUXTON_178      C       BUXTON  BUXTON_WUR_AseI_NsiI_final_run1 1       2       AseI    NsiI    3       3       Scabiosa columbaria
H53KHCCXY       5       CCTA    CCAG    WUR_178 C       WUR     BUXTON_WUR_AseI_NsiI_final_run1 2       2       AseI    NsiI    3       3       Scabiosa columbaria
H53KHCCXY       5       TTAC    CCAG    BUXTON_169      C       BUXTON  BUXTON_WUR_AseI_NsiI_final_run1 3       2       AseI    NsiI    3       3       Scabiosa columbaria
H53KHCCXY       5       AGGC    CCAG    WUR_169 C       WUR     BUXTON_WUR_AseI_NsiI_final_run1 4       2       AseI    NsiI    3       3       Scabiosa columbaria
H53KHCCXY       5       GAAGA   CCAG    BUXTON_175      SD      BUXTON  BUXTON_WUR_AseI_NsiI_final_run1 5       2       AseI    NsiI    3       3       Scabiosa columbaria
H53KHCCXY       5       CCTTC   CCAG    WUR_175 SD      WUR     BUXTON_WUR_AseI_NsiI_final_run1 6       2       AseI    NsiI    3       3       Scabiosa columbaria
```

## Start the pipeline

- Dry-run: 
	- `snakemake -n --use-conda`
- Run the pipeline: 
	`snakemake -j <threads> -p --use-conda` Replace <threads> by number of CPU???s to use on your server, e.g. 12

## Explanation of files in the output directory

It follows a description of all output files. Files that are important for downstream analysis are highlighted in bold. Files or Directories in italics are specific for the de-novo and reference branch respectively. 

-  __report.html__: A report summarizing all stats from the epiGBS analysis. Absolutely crucial to determine, whether your analysis run successfully or not.
- __multiQC_report.html__: A report summarizing fastQC files for 
- output_demultiplex: 
	- barcode_stacks.tsv: barcode file converted to the required stacks format  
	- clone: Directory containing read files that were PCR duplicates were removed  
	- clone-stacks: Directory containing demultiplexed and Watson-Crick separated read files. File names start with sample names. Files containing *rem* contain reads that failed the RAD-tag check.  
	- crick_R1.fq.gz: All filtered and demultiplexed crick forward reads, input for following steps  
	- crick_R2.fq.gz: All filtered and demultiplexed crick reverse reads, input for the following steps  
	- demultiplex.log: Log file containing stats from demultiplexing, also contains information about eventually discovered ???de novo??? barcodes  
	- Watson_R1.fq.gz: All filtered and demultiplexed watson forward reads, input for the following steps    
	- Watson_R2.fq.gz: All filtered and demultiplexed watson reverse reads, input for the following steps  
- *Output_denovo*: Containing all files created during de-novo reference creation. Skipped, if running in reference-mode.
	- Assembled.fq.gz: All read pairs (forward/R1 and reverse/R2), that were assembled/merged using Pear
	- Unassembled.R1.fq.gz: All forward/R1 reads that could not be assembled/merged  
	- Unassembled.R2.fq.gz: All forward/R2 reads that could not be assembled/merged  
	- consensus.fa: sequence file. outputfile from second clustering step, in which binary watson and crick reads are matched.    
	- consensus_cluster.fa: De novo reference sequence file. Outputfile from third clustering step, in which sequences from previous steps are clustered based on identity  
	- __consensus_cluster.renamed.fa__: sequences are identical to consensus_cluster.fa but fasta names are renamed. Input for mapping.
	- consensus_cluster.renamed.fa.fai: Indexed de novo reference sequences. 
	- make_reference.log: Log file containing some stats about the clustering statistics.  
- mapping, *output_mapping*: 
	- STAR_{joined,merged}_{crick,watson}: Directory, which contains the (indexed) reference converted to a three-letter alphabet. If working with a reference genome, differentiation between joined/merged absent. 
	- {crick, watson}.bam: Alignment file of crick or watson reads against the crick- or watson-converted reference, created by STAR
	- {crick,watson}.bam.bai: indexed alignment files, used as input for samtools        
	- {crick,watson}.vcf.gz: created by samtools, contains all variant positions in {crick,watson} reads against the reference for each sample          
	- merged.tsv.gz: created by a custom script that merges crick.vcf.gz and watson.vcf.gz   
	- __snp.vcf.gz__: created by a custom script that extracts all SNPs from merged.tsv.gz. You can read this file using zcat or after unzipping it using gunzip
	- snp.vcf.gz.tbi: samtools tabix indexed SNP file            
	- header.sam    
	- mapping_variantcalling.log: Contains the mapping statistics 
	- __methylation.bed__: created by a custom script that extracts all methylated positions for each sample
	- heatmap.igv: Input for IGV (https://software.broadinstitute.org/software/igv/)  to visualize methylated positions in the genome  
- *Trimmed*: Directory containing adapter trimmed reads and their fastqc output files. For running in reference-mode only.
- multiQC_report_data: Directory containing all log and intermediary files that are created by MultiQC
- Fastqc: Directory containing individual fastQC reports
- Log: Directory containing log-files

## When not to run the pipeline?

If you want to determine methylation in restriction enzyme overhang. The original sequence will be replaced by the expected during demultiplexing if a sequence error occurs. This does not take into account C-T conversions.

## Quality control or "How to discover errors?"

Recommendation: Run fastq-screen on raw data to determine sources of contamination (e.g. by sharing a lane with other customers, human DNA, phiX, vectors and adapters). https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html and https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_documentation.html

### MultiQC report:

- MultiQC wrappes the FastQC report of several samples and gives you the posibility to compare the quality plots of all or specific samples directly with each other. All functions of MultiQC are explained in a tutorial video that is provided as a link in the report.
- you can consider the FastQC documentation to interpretate most of the FastQC plots: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/. However, epiGBS libraries have the following specific characteristics that will be reflected in the quality reports. You can also compare this with the RRBS example from the FastQC website: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/RRBS_fastqc.html
	- The sequence duplication levels will be higher than average because restriction enzymes were used to create the libraries versus random shearing for e.g. whole genome sequencing.
	- The per base sequence content of the first nucleotides will re-construct the overhang of the used restriction enzyme after demultiplexing the reads. The ???C??? content is usually low and ???T??? is high.

### EpiGBS-report:

- Number of clones: Check the number and distributions of clones. The percentage of clone reads should be < XX%
- Demultiplexing: The number of reads per sample should be comparable and agree with the expected numbers. The amount of Watson and Crick reads for a specific sample should be similar, too. The optimal number of reads depends on the expected number of fragments (> 10 reads per fragment). 
- De novo reference: 
	- Check the number of assembled reads. Usually the percentage of assembled reads should be higher than unassembled reads. However, it willdepend from the size of the fragments that you expect. If the fragment size is greater than twice the chosen read site, the amount of assembled reads will be small.
	- The number of consenus clusters (third clustering step) should be comparable with the numbers of fragments that you expected, e.g. by performing an in silico digest using R packages like SimRAD{Updating}
- Mapping:
	- The mapping percentage for all reads should be higher than 40%. In general, we observe a lower mapping percentage for assembled than for joined reads.
- SNP calling
	- The amount of detected SNPs will depend from the used species and the diversity of samples. The SNP depth should be greater than 10 for reliable calling. 
- Methylation calling
	- The amount of methylated cytosines and their context will depend from the used species and used restriction enzymes. In general, the depth should be greater than 10 for reliable calling. 

## Fix errors 

### Clone percentage 

#### Problem:
I have a very high percentage of clone reads

#### Fix:

- Check the length of the Wobble in your adapter sequence and the set number in the barcode file. 
- Check the quality and quantity of input DNA, as well as the expected length of DNA fragments.

### Demultiplexing 

#### Problem:
One or more samples have small amounts of recovered reads.

#### Fix: 

- Check the labwork
- Barcode log file (output.dir/output-demultiplex/clone-stacks/process_radtags.clone.log) contains de-novo discovered barcodes. Complete set of Watson and Crick, reverse and forward? ??? check your barcode file and experimental design file again
- Check DNA quantification and quality of the missing samples
- RAD-tag: Pipeline checks for the presence of the restriction enzyme overhang (RAD-tag). If this sequence contains ???C???s???, they can be converted to ???T???s???. The RAD-tag check allows one nucleotide mismatch. Depending from the chosen enzyme combination, two or more mismatches should be allowed. Possible fix is to switch off RAD-tag check by opening demultiplex/barcode_stacks.py and change line XX from cmd += "-r -D --inline_inline " to cmd += "-r -D --inline_inline --disable_rad_check "
- Check the .rm. files in output.dir/output-demultiplex/clone-stacks/ in the QC. Do they show a common sequence in the beginning of the read different from your expected RAD-tag?

#### Problem:
The coverage is too low in the methylation bed file and after filtering on coverage (>10) only few positions remain.

#### Fix:
The required amount of reads will depend from the expected number of DNA fragments. You can calculate this for your (or a related) species by using R packages like SimRAD{Updating} (https://cran.r-project.org/package=SimRAD)

## Example Config Files

### Reference

```
#path to output directory
output_dir: "/mnt/nfs/bioinfdata/home/NIOO/fleurg/projects/epigbs-snakemake-GT/output"
#input directory where raw reads are
input_dir       : "/mnt/nfs/bioinfdata/home/NIOO/fleurg/projects/epigbs-snakemake-GT/data"
barcodes: "barcodes.tsv"
#barcode file(barcode file should be kept inside input directory) and enzymes will be included in barcode file
tmpdir        : "/tmp"
#Projectname to ensure there is no overlap between tmpfiles in the demultiplexing step  
#number of sequencing cycles
cycles        : 150
#Sequence Reads
Read1 : "NIOO_DWGT04032019_FKDL190727766-1a_HYF2KCCXY_L2_1.fq.gz"
Read2 : "NIOO_DWGT04032019_FKDL190727766-1a_HYF2KCCXY_L2_2.fq.gz"
#genome directory (leave it blank in denovo mode)
ref_dir: "/mnt/nfs/bioinfdata/home/NIOO/fleurg/projects/epigbs-snakemake-GT"
#genome name (leave it blank in denovo mode)
genome: "NC_031768.fasta"

#mode of running pipeline (set denovo or reference)
mode: "reference"
#set parameters or write default
#identity in decimal number e.g. for 90% identity write 0.90
param_denovo:
  identity: ""
  min-depth: ""
  max-depth: ""

  
#SNP calling parameters
param_SNPcalling:
  max-depth: ""
  min-MQ: ""
  min-BQ: ""

```

### De novo

```
#path to output directory
output_dir: "/mnt/nfs/bioinfdata/home/NIOO/fleurg/projects/epigbs-snakemake-denovo-190829/output"
#input directory where raw reads are
input_dir       : "/mnt/nfs/bioinfdata/home/NIOO/fleurg/epiGBS/data"
barcodes: "barcodes.tsv"
#barcode file(barcode file should be kept inside input directory) and enzymes will be included in barcode file

tmpdir        : "/tmp"
#Projectname to ensure there is no overlap between tmpfiles in the demultiplexing step
#number of sequencing cycles
cycles        : 150
#Sequence Reads
Read1 : "RRBS_KD17072296_H53KHCCXY_L5_1.fq.gz"
Read2 : "RRBS_KD17072296_H53KHCCXY_L5_2.fq.gz"
#genome directory (leave it  blank in denovo mode)
ref_dir: ""
#genome name (leave it blank in denovo mode)
genome: ""

#mode of running pipeline (set denovo or reference)
mode: "denovo"
#set parameters or write default
#identity in decimal number e.g. for 90% identity write 0.90
param_denovo:
  identity: ""
  min-depth: ""
  max-depth: ""


#SNP calling parameters
param_SNPcalling:
  max-depth: ""
  min-MQ: ""
  min-BQ: ""

```

### De novo with custom parameters

```
#path to output directory
output_dir: "/mnt/nfs/bioinfdata/home/NIOO/fleurg/projects/epigbs-snakemake-denovo-190829/output"
#input directory where raw reads are
input_dir       : "/mnt/nfs/bioinfdata/home/NIOO/fleurg/epiGBS/data"
barcodes: "barcodes.tsv"
#barcode file(barcode file should be kept inside input directory) and enzymes will be included in barcode file

tmpdir        : "/tmp"
#Projectname to ensure there is no overlap between tmpfiles in the demultiplexing step
#number of sequencing cycles
cycles        : 150
#Sequence Reads
Read1 : "RRBS_KD17072296_H53KHCCXY_L5_1.fq.gz"
Read2 : "RRBS_KD17072296_H53KHCCXY_L5_2.fq.gz"
#genome directory(leave it  blank in denovo mode)
ref_dir: ""
#genome name (leave it blank in denovo mode)
genome: ""

#mode of running pipeline (set denovo or reference)
mode: "denovo"
#set parameters or write default
#identity in decimal number e.g. for 90% identity write 0.90
param_denovo:
  identity: "0.85"
  min-depth: ""
  max-depth: ""


#SNP calling parameters
param_SNPcalling:
  max-depth: ""
  min-MQ: ""
  min-BQ: ""

```

## List of used software and references

### Software 

- [Snakemake 5.4.5](https://snakemake.readthedocs.io/en/stable/)
- [Conda](https://docs.conda.io/en/latest/index.html)
- [Stacks](http://catchenlab.life.illinois.edu/stacks/)
- [Python 3.7](https://www.python.org/)
- [R + R package to render Rmd](https://www.r-project.org/)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- [Trim-galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore)
- Pear
- [Seqtk](https://github.com/lh3/seqtk)
- STAR
- vsearch
- [Samtools](http://www.htslib.org/)

### References

1.	K??ster, J. & Rahmann, S. Snakemake???a scalable bioinformatics workflow engine. Bioinformatics 28, 2520???2522 (2012).
2.	Gr??ning, B. et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat. Methods 15, 475???476 (2018).
3.	Catchen, J., Hohenlohe, P. A., Bassham, S., Amores, A. & Cresko, W. A. Stacks: an analysis tool set for population genomics. Mol. Ecol. 22, 3124???3140 (2013).
4.	Catchen, J. M., Amores, A., Hohenlohe, P., Cresko, W. & Postlethwait, J. H. Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. G3 Genes Genomes Genet. 1, 171???182 (2011).
5.	Stacks 2: Analytical Methods for Paired-end Sequencing Improve RADseq-based Population Genomics | bioRxiv. Available at: https://www.biorxiv.org/content/10.1101/615385v1. (Accessed: 27th August 2019)
6.	Andrews, Simon. FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc. (2010).
7.	Zhang, J., Kobert, K., Flouri, T. & Stamatakis, A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30, 614???620 (2014).
8.	Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15???21 (2013).
9.	Rognes, T., Flouri, T., Nichols, B., Quince, C. & Mah??, F. VSEARCH: a versatile open source tool for metagenomics. PeerJ 4, (2016).
10.	Lepais, O. & Weir, J. T. SimRAD: an R package for simulation-based prediction of the number of loci expected in RADseq and similar genotyping by sequencing approaches. Mol. Ecol. Resour. 14, 1314???1321 (2014).
