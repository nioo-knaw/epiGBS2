
param_oligo=config["param_demultiplex"]["nOligo"]

#We always use 3 as the number of oligo for our adapters. (If we would ever change this we can account for this)
def getParam_oligo(param_oligo):
    if param_oligo == "default" or param_oligo == "":
        id = 3
    else:
        id = param_oligo
    return id

#Clone_filter removes PCR duplicates
#It compares reads that are removes reads that are completely identical 
#including the UMI/wobble/oligo's (random bases that are added at the beginning of the reads)
#as if all of these are identical it will most likely be a PCR artifact
rule clone_filter:
    input:
        barcodes=expand("{path}/{bar}", path=config["input_dir"], bar=config["barcodes"]),
        R1=expand("{path}/{{run}}_R1.fq.gz",path=config["input_dir"],R1=RUN),
        R2=expand("{path}/{{run}}_R2.fq.gz",path=config["input_dir"],R2=RUN)
    params:
        output_dir=expand("{path}/output_demultiplex/", path=config["output_dir"]),
        param_oligo=getParam_oligo(param_oligo)
    output:
        R1=expand("{path}/output_demultiplex/clone/{{run}}_R1.1.fq.gz",path=config["output_dir"]),
        R2=expand("{path}/output_demultiplex/clone/{{run}}_R2.2.fq.gz",path=config["output_dir"])
    conda:
        "../env/stacks.yaml"
    threads: 1
    shell: 
        "clone_filter -1 {input.R1} -2 {input.R2} -o {params.output_dir}/clone/ --oligo_len_1 {params.param_oligo} --oligo_len_2 {params.param_oligo} --inline_inline -i gzfastq"

#Stacks and the rest of the pipeline need to have specific files for barcodes, 
#the format is different and we add the control nucleotide and we need to split it per run so we can demultiplex them in parallel
#popmap (to which population do samples belong),
#and an input file for SNPFilter report to add nice colours to the plots based on a priori clustering
rule make_stacks_files:
    input:
        barcodes=expand("{path}/{bar}", path=config["input_dir"], bar=config["barcodes"])
    output:
        barcodes=expand("{path}/output_demultiplex/barcode_stacks{run}.tsv", path=config["output_dir"], bar=config["barcodes"],run=RUN)
    params:
        output_dir=expand("{path}/output_demultiplex/",path=config["output_dir"])
    conda:
        "../env/R.yaml"
    shell:
        "Rscript src/demultiplex/createFilesFrombarcodes.R {input.barcodes} {params.output_dir}"

#We then do some snakemake magic to run process radtags for each run.
#This outputs 2 files per sample (the R1 and R2) 
#TODO figure out to have this not be completely dependent on the process_radtags.log file being generated.
#Howver this is complicated due to the samples needing to be split up per file...
#TODO make this less janky... so we do not do weird shit with the logs...
#Figure out what happesn if this breaks.
rule process_radtags:
    input:
        barcodes=expand("{path}/output_demultiplex/barcodeStacks{{run}}.tsv", path=config["output_dir"], bar=config["barcodes"]),
        R1=expand("{path}/output_demultiplex/clone/{{run}}_R1.1.fq.gz",path=config["output_dir"]),
        R2=expand("{path}/output_demultiplex/clone/{{run}}_R2.2.fq.gz",path=config["output_dir"])
    params:
        output_dir=expand("{path}/output_demultiplex/logs/{{run}}/",path=config["output_dir"])
    conda:
        "../env/stacks.yaml"
    threads: THREADSPERRUN//1
    output:
        log=expand("{path}/output_demultiplex/logs/{{run}}/process_radtags.clone_filter.log",path=config["output_dir"])
    shell:
        "process_radtags -1 {input.R1} -2 {input.R2} -o {params.output_dir} -b {input.barcodes} --renz_1 aseI --renz_2 nsiI -c --inline-inline --threads {threads}"

#This moves all samples into the demultiplex/samples directory
rule moveDemultiplexFiles:
    input:
        log=expand("{path}/output_demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["output_dir"],run=RUN)
    params:
        input_dir=expand("{path}/output_demultiplex/logs/",path=config["output_dir"]),
        output_dir=expand("{path}/output_demultiplex/clone-stacks/",path=config["output_dir"]),
        log=expand("{path}/logs/",path=config["output_dir"])
    output:
        samplesR1=expand("{path}/output_demultiplex/clone-stacks/{samples}.1.fq.gz",path=config["output_dir"],samples=SAMPLES),
        samplesR2=expand("{path}/output_demultiplex/clone-stacks/{samples}.2.fq.gz",path=config["output_dir"],samples=SAMPLES)
    shell:
        """
        mv {params.input_dir}/*/*.fq.gz {params.output_dir}/
        """