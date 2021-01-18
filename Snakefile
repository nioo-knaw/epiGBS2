
configfile: "config.yaml"
import pandas as pd
import os
import random
df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
flowCell = df.Flowcell[0]
lane = df.Lane[0]
projectName=random.randint(1,1000000) #To ensure non overlapping tmp directories

if config["mode"] == "reference":
    include: "src/rules/fastqc-ref.rules"
    include: "src/rules/with_reference.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/report-ref.rules"

if config["mode"] == "denovo":
    include: "src/rules/denovo_reference.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/bismarkDenovo.rules"
    include: "src/rules/trimming.rules"
    include: "src/rules/snp_calling.rules"

if config["mode"]== "reference":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/multiQC_report.html \
            {out}/fastqc/ \
            {out}/trimmed/Watson_R1_trim_fastqc.html \
            {out}/trimmed/Crick_R2_trim_fastqc.html \
            {out}/report.html \
            {out}/mapping/watson.bam \
            {out}/mapping/crick.bam \
            {out}/mapping/methylation.bed \
            {out}/mapping/snp.vcf.gz".split(),out=config["output_dir"])

if config["mode"]== "denovo":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/output_denovo/consensus_cluster.renamed.fa \
		    {out}/log/{sample}_read-info.txt \
		    {out}/log/{sample}_untrimmed_filt_read-info.txt \
		    {out}/log/{sample}_trimmed_three_read-info.txt \
		    {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
            {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
		    {out}/alignment/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.bam \
		    {out}/methylation_calling/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.CX_report.txt \
		    {out}/methylation_calling/{sample}_bismark.cov \
            {out}/snp_calling/{sample}_snp.raw.vcf.out",path=config["output_dir"]).split(),out=config["output_dir"],sample=SAMPLES)


# if config["mode"]== "denovo":
#     rule all:
#         input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
#             {out}/output_demultiplex/Watson_R1.fq.gz \
#             {out}/output_demultiplex/Crick_R2.fq.gz \
#             {out}/output_demultiplex/Crick_R1.fq.gz \
#             {out}/fastqc/ \
#             {out}/multiQC_report.html \
#             {out}/report.html \
#             {out}/output_denovo/consensus_cluster.renamed.fa \
# 		    {out}/log/{sample}_read-info.txt \
# 		    {out}/log/{sample}_untrimmed_filt_read-info.txt \
# 		    {out}/log/{sample}_trimmed_three_read-info.txt \
# 		    {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
#             {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
# 		    {out}/alignment/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.bam \
# 		    {out}/methylation_calling/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.CX_report.txt \
# 		    {out}/methylation_calling/{sample}_bismark.cov".split(), out=config["output_dir"],sample=SAMPLES)
# if config["mode"] == "denovo":
#     include: "src/rules/denovo_reference.rules"
#     include: "src/rules/demultiplex.rules"
#     include: "src/rules/fastqc.rules"
#     include: "src/rules/report.rules"
#     include: "src/rules/bismarkDenovo.rules"
#     include: "src/rules/trimming.rules"
