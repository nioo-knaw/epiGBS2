
configfile: "config.yaml"


if config["mode"] == "reference":
    include: "src/rules/fastqc-ref.rules"
    include: "src/rules/with_reference.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/report-ref.rules"

if config["mode"] == "denovo":
    include: "src/rules/denovo_reference.rules"
    include: "src/rules/demultiplex.rules"
    include: "src/rules/fastqc.rules"
    include: "src/rules/report.rules"

if config["mode"]== "reference":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/multiQC_report.html \
            {out}/fastqc/ \
            {out}/trimmed/Watson_R1_val_1_fastqc.html \
            {out}/trimmed/Crick_R2_val_2_fastqc.html \
            {out}/report.html \
            {out}/mapping/watson.bam \
            {out}/mapping/crick.bam \
            {out}/mapping/methylation.bed \
            {out}/mapping/snp.vcf.gz".split(), out=config["output_dir"])


if config["mode"]== "denovo":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/fastqc/ \
            {out}/multiQC_report.html \
            {out}/report.html \
            {out}/mapping/watson.bam \
            {out}/mapping/crick.bam \
            {out}/output_denovo/consensus_cluster.renamed.fa \
            {out}/mapping/methylation.bed \
            {out}/mapping/snp.vcf.gz".split(), out=config["output_dir"])
