
configfile: "config.yaml"


if config["mode"] == "reference":
    include: "rules/fastqc-ref.rules"
    include: "rules/with_reference.rules"
    include: "rules/demultiplex.rules"
    include: "rules/report-ref.rules"

if config["mode"] == "denovo":
    include: "rules/denovo_reference.rules"
    include: "rules/demultiplex.rules"
    include: "rules/fastqc.rules"
    include: "rules/report.rules"

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
            {out}/output_mapping/methylation.bed \
            {out}/output_mapping/snp.vcf.gz".split(), out=config["output_dir"])


if config["mode"]== "denovo":
    rule all:
        input: expand("{out}/output_demultiplex/Watson_R2.fq.gz \
            {out}/output_demultiplex/Watson_R1.fq.gz \
            {out}/output_demultiplex/Crick_R2.fq.gz \
            {out}/output_demultiplex/Crick_R1.fq.gz \
            {out}/fastqc/ \
            {out}/multiQC_report.html \
            {out}/report.html \
            {out}/output_denovo/consensus_cluster.renamed.fa \
            {out}/mapping/methylation.bed \
            {out}/mapping/snp.vcf.gz".split(), out=config["output_dir"])
