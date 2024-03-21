
configfile: "config.yaml"
from snakemake.utils import Paramspace
import pandas as pd
import os
import random
df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
SAMPLE = SAMPLES[0]
flowCell = "Redudant"
lane = "Redundant"
projectName=random.randint(1,1000000) #To ensure non overlapping tmp directories
RAWREADSR1 = df.rawR1.str.replace(".fq.gz","",regex=False).unique()
RAWREADSR2 = df.rawR2.str.replace(".fq.gz","",regex=False).unique()
RUN = df.rawR1.str.replace("_R1.fq.gz","",regex=False).unique()
OLIGOR1 = df.Wobble_R1.unique()
OLIGOR1 = df.Wobble_R2.unique()

paramspace = Paramspace(pd.read_csv("src/parameter_test/paramTest.tsv", sep="\t"))



if config["mode"] == "reference":
    include: "src/rules/fastqc-ref.rules"
    include: "src/rules/reference.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/trimming.rules"


if config["mode"] == "denovo":
    include: "src/rules/demultiplex.rules"
    include: "src/rules/denovo.rules"
    include: "src/rules/trimming.rules"
    include: "src/rules/fastqc.rules"

if config["mode"] == "paramTest":
    include: "src/rules/demultiplex.rules"
    include: "src/rules/paramTest.rules"
    include: "src/rules/trimming.rules"


if config["mode"]== "reference":
    rule all:
        input: expand("{out}/output_demultiplex/clone-stacks/{sample}-Watson.1.fq.gz \
            {out}/output_demultiplex/clone-stacks/{sample}-Watson.2.fq.gz \
            {out}/output_demultiplex/clone-stacks/{sample}-Crick.1.fq.gz \
            {out}/output_demultiplex/clone-stacks/{sample}-Crick.2.fq.gz \
            {out}/fastqc/ \
            {out}/multiQC_report.html \
		    {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
            {out}/cutadapt/{sample}_trimmed_filt_merged.2.fq.gz \
		    {out}/alignment/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.bam \
		    {out}/methylation_calling/{sample}_bismark_bt2_pe.CX_report.txt.gz \
            {out}/snp_calling/snp.vcf.gz".split(),out=config["output_dir"],sample=SAMPLES)

if config["mode"]== "denovo":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/fastqc/ \
            {out}/multiQC_report.html \
            {out}/output_denovo/consensus_cluster.renamed.fa \
		    {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
            {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
		    {out}/alignment/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.bam \
		    {out}/methylation_calling/{sample}_bismark_bt2_pe.CX_report.txt.gz \
            {out}/snp_calling/snp.vcf.gz".split(),out=config["output_dir"],sample=SAMPLES)

if config["mode"]== "paramTest":
    rule all:
        input: expand("{out}/output_demultiplex/clone-stacks/{samples}-Watson.1.fq.gz \
            {out}/output_demultiplex/clone-stacks/{samples}-Watson.2.fq.gz \
            {out}/output_demultiplex/clone-stacks/{samples}-Crick.1.fq.gz \
            {out}/output_demultiplex/clone-stacks/{samples}-Crick.2.fq.gz \
            {out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/paramTest/{params}/alignment/{sample}_trimmed_filt_merged.1_bismark_bt2_pe.bam \
            {out}/paramTest/{params}/output_denovo/consensus_cluster.renamed.fa \
            {out}/paramTest/Assembled.txt \
            {out}/paramTest/averageDepth.txt \
            {out}/paramTest/denovoParameter.tsv \
            {out}/paramTest/denovoParameter.tiff \
		    {out}/cutadapt/{sample}_trimmed_filt_merged.1.fq.gz \
            {out}/cutadapt/{sample}_trimmed_filt_merged.2.fq.gz".split(),out=config["output_dir"],params=paramspace.instance_patterns,sample=SAMPLE,samples=SAMPLES)
    

if config["mode"] == "legacy":
    include: "src/rules/legacy.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/report.rules"


if config["mode"]== "legacy":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/mapping/watson.bam \
            {out}/mapping/crick.bam \
            {out}/output_denovo/consensus_cluster.renamed.fa \
            {out}/mapping/methylation.bed \
            {out}/mapping/snp.vcf.gz".split(), out=config["output_dir"])


