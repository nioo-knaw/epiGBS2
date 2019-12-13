#!/usr/bin/env python
__author__ = 'thomasvangurp'
# Date created: 22/11/2014 (europe date)
# Function: Pipeline for mapping reads to reference
#Python version: 3.6.3
#External dependencies: samtools,pysam,methylation_calling.py SNP_calling.py
#Known bugs: None
#Modifications: None
import pysam
import argparse
import subprocess
import tempfile
import os
import shutil
import sys
import time
from Bio import SeqIO
from Bio import Restriction
from os import system


def getScriptPath():
    return os.path.dirname(__file__)

def parse_args():
    "Pass command line arguments"
    if not sys.argv[1:]:
        sys.argv.append('-h')
    parser = argparse.ArgumentParser(description='use bwameth for mapping reads')
    #input files
    parser.add_argument('-s','--sequences',
                        help='number of sequences to take for testing')
    parser.add_argument('--max_depth',default=999999999,
                        help='maximum depth for SNP calling')
    parser.add_argument('--min_MQ',default=0,
                        help='minimum mapQ')
    parser.add_argument('--min_BQ',default=15,
                        help='minimum baseQ/BAQ')
    parser.add_argument('--subsample_treshold',
                        help='Subsample treshold',default='100000')
    parser.add_argument('--tmpdir',
                        help='tmp directory',default="/tmp/")
    parser.add_argument('--input_dir',
                        help='optional: Choose input directory')
    parser.add_argument('--reads_R1',
                    help='Forward unmerged reads')
    parser.add_argument('--reads_R2',
                        help='Reverse unmerged reads')
    parser.add_argument('--merged',
                        help='merged watson and crick fastq')
    parser.add_argument('--reference',
                    help='reference clusters')
    parser.add_argument('--refgenome',
                        help='reference genome instead of clusters')
    parser.add_argument('-b','--barcodes',
                    help='Barcodes used in output')
    parser.add_argument('--species',
                        help='Species: if selected only that species will be putin BAM RG header')
    parser.add_argument('--bamout',
                        help='output for bam file with RGs')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--log',
                        help='log of output operation')
    parser.add_argument('--output_dir',
                        help='optional: Choose output directory')
    parser.add_argument('--watson_vcf',
                        help='watson vcf output')
    parser.add_argument('--crick_vcf',
                        help='crick vcf output')
    parser.add_argument('--snp_vcf',
                        help='vcf output snp')
    parser.add_argument('--methylation_vcf',
                        help='Methylation vcf output')
    parser.add_argument('--heatmap',
                        help='heatmap output methylation')
    args = parser.parse_args()
    if args.input_dir:
        if not args.reads_R1:
            args.reads_R1 = os.path.join(args.input_dir,'Unassembled.R1.watson.fq.gz')
        if not args.reads_R2:
            args.reads_R2 = os.path.join(args.input_dir,'Unassembled.R2.crick.fq.gz')
        if not args.merged:
            args.merged = os.path.join(args.input_dir,'Assembled.fq.gz')
        if args.reference == None and args.refgenome == None:
            args.reference = os.path.join(args.input_dir,'consensus_cluster.renamed.fa')
        if args.barcodes == None:
            args.barcodes = os.path.join(args.input_dir,'barcodes.csv')
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        if not args.log:
            args.log = os.path.join(args.output_dir,'mapping_variantcalling.log')
        args.watson_vcf = os.path.join(args.output_dir,'watson.vcf.gz')
        args.crick_vcf = os.path.join(args.output_dir,'crick.vcf.gz')
        args.snp_vcf = os.path.join(args.output_dir,'snp.vcf.gz')
        args.methylation_vcf = os.path.join(args.output_dir,'methylation.vcf.gz')
        args.heatmap = os.path.join(args.output_dir,'heatmap.igv')
        args.mastermeth = os.path.join(args.output_dir,'methylation.bed')
    return args

def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log,'a') as log:
        log.write("now starting:\t%s\n"%log_message)
        log.write('running:\t%s\n'%(' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode().replace('\r','\n')
        stderr = stderr.decode().replace('\r','\n')
        if stdout:
            log.write('stdout:\n%s\n'%stdout)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        return_code = p.poll()
        if return_code:
            raise RuntimeError(stderr)
        log.write('finished:\t%s\n\n'%log_message)
    return 0


def run_STAR(in_files, args):
    "run STAR for mapping"

    in_files['bam_out'] = {}
    in_files['bam_out']['watson'] = os.path.join(args.output_dir, 'watson.bam')
    in_files['bam_out']['crick'] = os.path.join(args.output_dir, 'crick.bam')
    in_files['header'] = os.path.join(args.output_dir, 'header.sam')
    cmd = ["python src/mapping_varcall/map_STAR.py",
           '--reads_R1 %s' % args.reads_R1,
           '--reads_R2 %s' % args.reads_R2,
           '--merged %s' % args.merged,
           "--barcodes %s" % args.barcodes,
           "--tmpdir %s" % args.tmpdir,
           "--threads %s" % args.threads,
           "--output_dir %s" % args.output_dir]
    if not args.reference:
        cmd += ['--refgenome %s' % args.refgenome]
    else:
        cmd += ['--reference %s' % args.reference]
    if args.sequences != None:
        cmd += ['--sequences %s' % args.sequences]
    log = "Map reads using STAR"
    run_subprocess(cmd, args, log)

    return in_files


def variant_calling_samtools(in_files,args):
    """Do variant calling with freebayes"""
    #run mpileup on watson bam file
    in_files['vcf_out'] = {}
    in_files['vcf_out']['watson'] = os.path.join(args.output_dir,'watson.vcf.gz')
    in_files['vcf_out']['crick'] = os.path.join(args.output_dir,'crick.vcf.gz')

    cmd = ["samtools mpileup --reference %s -gt DP,AD,INFO/AD" % (args.reference) +
           " --max-depth  %s " %(args.max_depth) +  # call at max-depth of 10.000.000
           "-q %s " %(args.min_MQ) +  # Do not skip alignments with low mapQ
           "-Q %s " %(args.min_BQ) +  # Skip bases with baseQ/BAQ smaller than 15
           "--skip-indels " +  # skip indels
           "-vu %s" % (
           in_files['bam_out']['watson']) +  # v = generate genotype likelihoods in VCF format u = uncompressed
           "|grep -v '^##contig='|bgzip -c > %s" % (in_files['vcf_out']['watson'])]

    log = "use samtools mpileup to get variant observation counts for watson"
    run_subprocess(cmd, args, log)


    cmd = ["samtools mpileup --reference %s -gt DP,AD,INFO/AD" % (args.reference) +
           " --max-depth  %s " %(args.max_depth) + #call at max-depth of 10.000.000
           "-q %s " %(args.min_MQ) + #Do not skip alignments with low mapQ #TODO: investigate option
           "-Q %s " %(args.min_BQ) + #Skip bases with baseQ/BAQ smaller than 15
           "--skip-indels " + #skip indels
           "-vu %s" % (in_files['bam_out']['crick']) + #v = generate genotype likelihoods in VCF format u = uncompressed
           "|grep -v '^##contig='|bgzip -c > %s" % (in_files['vcf_out']['crick'])]

    log = "use samtools mpileup to get variant observation counts for crick"
    run_subprocess(cmd, args, log)
    return in_files

def merge_watson_crick(in_files, args):
    """create merged.tsv.gz with watson and crick calls merged"""
    if 'vcf_out' not in in_files:
        in_files['vcf_out'] = {}
        in_files['vcf_out']['watson'] = os.path.join(args.output_dir, 'watson.vcf.gz')
        in_files['vcf_out']['crick'] = os.path.join(args.output_dir, 'crick.vcf.gz')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir,'merged.tsv')
    cmd = ["python src/mapping_varcall/merge_watson_crick.py",
           "-w %s" % in_files['vcf_out']['watson'],
           "-c %s" % in_files['vcf_out']['crick'],
           "-o %s" % in_files['vcf_out']['merged']]

    log = "Create custom tsv file for combining watson and crick observation counts per individual"
    run_subprocess(cmd, args, log)
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    return in_files

def SNP_calling(in_files, args):
    """run SNP calling"""
    if 'vcf_out' not in in_files:
        in_files['vcf_out'] = {}
    in_files['vcf_out']['SNP'] = os.path.join(args.output_dir, 'snp.vcf')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    cmd = ["python src/mapping_varcall/SNP_calling.py",
           "-m %s" % in_files['vcf_out']['merged'],
           "-s %s" % in_files['vcf_out']['SNP'],
           "-w %s" % os.path.join(args.output_dir, 'watson.vcf.gz')]
    log = "perform SNP calling"
    run_subprocess(cmd, args, log)

    return in_files


def methylation_calling(in_files,args):
    "run methylation calling script."
    log = ["Run methylation calling script"]
    in_files['vcf_out']['SNP'] = os.path.join(args.output_dir, 'snp.vcf.gz')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    cmd = ["python src/mapping_varcall/methylation_calling.py",
           " -r %s"%(args.reference),
           " -m %s"%(in_files['vcf_out']['merged']),
           " -s %s"%(in_files['vcf_out']['SNP']),
           " -o %s"%(os.path.join(args.output_dir,'methylation.bed')),
           " -heat %s"%(os.path.join(args.output_dir,'heatmap.igv'))
           ]



    run_subprocess(cmd,args,log)
    return in_files

def main():
    "Main function loop"
    args = parse_args()
    #Make sure log is empty at start
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: discover files in input #todo
    files = {}
    #Step 2: map reads using STAR
    #TODO: replace for running map_STAR
    files = run_STAR(files,args)
    if args.refgenome:
        args.reference = args.refgenome
    files = variant_calling_samtools(files, args)
    files = merge_watson_crick(files,args)
    files = SNP_calling(files, args)
    files = methylation_calling(files,args)
    print('done')
if __name__ == '__main__':
    main()
