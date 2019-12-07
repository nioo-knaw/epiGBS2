#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
import argparse
import subprocess
import os
import math
import gzip
import tempfile
import urllib
import json
import ssl

__author__ = 'thomasvangurp'
__description__ = "map reads orders of magnitudes faster using STAR"


def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='use STAR for mapping reads')
    #input files
    parser.add_argument('-s', '--sequences',
                        help='number of sequences to take for testing, subsamples the input read files')
    parser.add_argument('--tmpdir',
                        help='tmp directory',default="/tmp/")
    parser.add_argument('--input_dir',
                        help='optional: Choose input directory with trimmed reads')
    parser.add_argument('--crick_val_r1',
                        help='Crick trimmed forward reads')
    parser.add_argument('--crick_val_r2',
                        help='Crick trimmed reverse reads')
    parser.add_argument('--watson_val_r1',
                        help='Watson trimmed forward reads')
    parser.add_argument('--watson_val_r2',
                        help='Watson trimmed reverse reads')
    parser.add_argument('--refgenome',
                        help='reference genome')
    parser.add_argument('--barcodes',
                    help='Barcodes used in output')
    parser.add_argument('--species',
                        help='Species: if selected only that species will be put in BAM RG header')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--output_dir',
                        help='Choose output directory')
    parser.add_argument('--extraflags',
                        help='extra flags for testing')
    args = parser.parse_args()
    if args.input_dir:
        if not args.crick_val_r1:
            args.crick_val_r1 = os.path.join(args.input_dir,'Crick_R1_val_1.fq.gz')
        if not args.crick_val_r2:
            args.crick_val_r2 = os.path.join(args.input_dir,'Crick_R2_val_2.fq.gz')
        if not args.watson_val_r1:
            args.watson_val_r1 = os.path.join(args.input_dir,'Watson_R1_val_1.fq.gz')
        if not args.watson_val_r2:
            args.watson_val_r2 = os.path.join(args.input_dir,'Watson_R2_val_2.fq.gz')
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        if 'log' not in args:
            args.log = os.path.join(args.output_dir,'mapping_variantcalling.log')
        args.watson_vcf = os.path.join(args.output_dir,'watson.vcf')
        args.crick_vcf = os.path.join(args.output_dir,'crick.vcf')
        args.snp_vcf = os.path.join(args.output_dir,'snp.vcf')
        args.methylation_vcf = os.path.join(args.output_dir,'methylation.vcf')
        args.heatmap = os.path.join(args.output_dir,'heatmap.igv')
        #2 bed files should be made for subsequent analysis using Rnbeads or other software
        args.mastermeth = os.path.join(args.output_dir,'methylation.bed')
    args.tmpdir = tempfile.mkdtemp(suffix='STAR', prefix='tmp', dir=args.tmpdir)
    return args

def get_version():
    """get version of current script"""
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    while True:
        if '.git' in os.listdir(parent_dir):
            break
        parent_dir = os.path.dirname(parent_dir)
    git_log = os.path.join(parent_dir,'.git','logs','HEAD')
    handle = open(git_log,'r')
    log_lines = [l.split('\t') for l in handle.readlines()]
    #now get latest github commit
    url = 'https://api.github.com/repos/thomasvangurp/epiGBS/commits'
    context = ssl._create_unverified_context()
    result = json.load(urllib.urlopen(url,context=context))
    print('')



def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log, 'a') as log:
        log.write("now starting:\t%s\n" % log_message)
        log.write('running:\t%s\n' % (' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode().replace('\r', '\n')
        stderr = stderr.decode().replace('\r', '\n')
        if stdout:
            log.write('stdout:\n%s\n' % stdout)
        if stderr:
            log.write('stderr:\n%s\n' % stderr)
        return_code = p.poll()
        if return_code:
            raise RuntimeError(stderr)
        log.write('finished:\t%s\n\n' % log_message)
    return 0



def process_reads_watson(args):
    """process watson trimmed reads and make them ready for mapping with STAR"""
    watson_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_r1', dir=args.tmpdir,
                                                   delete=False)
    watson_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_r2', dir=args.tmpdir,
                                                   delete=False)
    args.watson_r1 = watson_r1.name
    args.watson_r2 = watson_r2.name
    print('Started processing watson reads')
    if args.watson_val_r1.endswith('.gz'):
        w_r1_handle = gzip.open(args.watson_val_r1, 'rt')
        w_r2_handle = gzip.open(args.watson_val_r2, 'rt')
    else:
        w_r1_handle = open(args.watson_val_r1, 'rt')
        w_r2_handle = open(args.watson_val_r2, 'rt')
    #make 4 file handles for forward and reverse watson and crick
    watson_r1_handle = open(args.watson_r1, 'w')
    watson_r2_handle = open(args.watson_r2, 'w')
    j = 0
    while True:
        w_r1 = []
        w_r2 = []
        for i in range(4):
            try:
                w_r1.append(next(w_r1_handle))
                w_r2.append(next(w_r2_handle))
            except StopIteration:
                break
        j += 1
        try:
            if int(args.sequences) == j:
                break
        except TypeError:
            pass
        if not j % 1000000:
            print('Processed %s reads' % (j))
        if not w_r1:
            break
        convert_w_r1 = w_r1[1].upper().replace('C', 'T')
        convert_w_r2 = w_r2[1].upper().replace('G', 'A')
        c_pos_w = [str(n) for n, i in enumerate(w_r1[1]) if i.upper() == 'C']
        g_pos_w = [str(n) for n, i in enumerate(w_r2[1].rstrip()[::-1]) if i.upper() == 'G']
        header_w = '@%s' % (w_r1[0][1:-1].replace(' ', '|').replace('\t', '|'))
        header_w += '|%s\n' % (','.join(c_pos_w) + '|' + ','.join(g_pos_w))
        watson_r1_handle.write(header_w + convert_w_r1 + '+\n' + w_r1[3])
        #print(read_r1[3])
        watson_r2_handle.write(header_w + convert_w_r2 + '+\n' + w_r2[3])
    watson_r1_handle.close()
    watson_r2_handle.close()
    return args

def process_reads_crick(args):
    """process crick trimmed reads and make them ready for mapping with STAR"""

    crick_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_r1', dir=args.tmpdir,
                                                   delete=False)
    crick_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_r2', dir=args.tmpdir,
                                                   delete=False)
    args.crick_r1 = crick_r1.name
    args.crick_r2 = crick_r2.name

    print('Started processing crick reads')
    if args.crick_val_r1.endswith('.gz'):
        c_r1_handle = gzip.open(args.crick_val_r1, 'rt')
        c_r2_handle = gzip.open(args.crick_val_r2, 'rt')
    else:
        c_r1_handle = open(args.crick_val_r1, 'rt')
        c_r2_handle = open(args.crick_val_r2, 'rt')
    #make 4 file handles for forward and reverse watson and crick
    crick_r1_handle = open(args.crick_r1, 'w')
    crick_r2_handle = open(args.crick_r2, 'w')
    j = 0
    while True:
        c_r1 = []
        c_r2 = []
        for i in range(4):
            try:
                read_r1.append(next(r1_handle))
                read_r2.append(next(r2_handle))
            except StopIteration:
                break
        j += 1
        try:
            if int(args.sequences) == j:
                break
        except TypeError:
            pass
        if not j % 1000000:
            print('Processed %s reads' % (j))
        if not c_r1:
            break
        convert_c_r1 = c_r1[1].upper().replace('G', 'A')
        convert_c_r2 = c_r2[1].upper().replace('C', 'T')
        g_pos_c = [str(n) for n, i in enumerate(c_r1[1]) if i.upper() == 'G']
        c_pos_c = [str(n) for n, i in enumerate(c_r2[1].rstrip()[::-1]) if i.upper() == 'C']
        header_c = '@%s' % (c_r1[0][1:-1].replace(' ', '|').replace('\t', '|'))
        header_c += '|%s\n' % (','.join(c_pos_c) + '|' + ','.join(g_pos_c))
        crick_r1_handle.write(header_c + convert_c_r1 + '+\n' + c_r1[3])
        #print(read_r1[3])
        crick_r2_handle.write(header_c + convert_c_r2 + '+\n' + c_r2[3])
    crick_r1_handle.close()
    crick_r2_handle.close()
    return args


def index_STAR(args):
    """make STAR index for reference genome"""
       # make STAR index folder for merged path
    ref_STAR_watson_index = os.path.join(args.output_dir,'STAR_watson')
    ref_STAR_crick_index = os.path.join(args.output_dir,'STAR_crick')
    if not os.path.exists(ref_STAR_watson_index):
        os.mkdir(ref_STAR_watson_index)
        os.mkdir(ref_STAR_crick_index)
    #ref_watson = os.path.join(ref_STAR_watson_index, '%s.ref.watson.fa' % args.species)
    #ref_crick = os.path.join(ref_STAR_crick_index, '%s.ref.crick.fa' % args.species)
    ref_watson = os.path.join(ref_STAR_watson_index,'ref.watson.fa')
    ref_crick = os.path.join(ref_STAR_crick_index,'ref.crick.fa')
    file_handle = open(args.refgenome, 'r')
    try:
        file_handle = open(args.refgenome, 'r')
        ref_len = 0
        ref_count = 0
        ref_watson_handle = open(ref_watson, 'w')
        ref_crick_handle = open(ref_crick, 'w')
        seq = ''
        for line in file_handle:
            if line.startswith('>'):
                if seq != '':
                    ref_len += len(seq)
                    ref_count += 1
                    ref_watson_handle.write(header + seq.upper().replace('C', 'T') + '\n')
                    ref_crick_handle.write(header + seq.upper().replace('G', 'A') + '\n')
                seq = ''
                header = line
            else:
                seq += line.rstrip('\n')
        ref_len += len(seq)
        ref_count += 1
        ref_watson_handle.write(header + seq.upper().replace('C', 'T') + '\n')
        ref_crick_handle.write(header + seq.upper().replace('G', 'A') + '\n')
    except Message:
        print("no reference genome found")
    # close file handles
    ref_watson_handle.close()
    ref_crick_handle.close()

    #iterate over input lines and write to references

    #MAKE LIST for indexes to be made
    index_list = [(ref_len, ref_count, ref_STAR_watson_index, ref_watson),
                  (ref_len, ref_count, ref_STAR_crick_index, ref_crick)]
    #calculate parameters for indexing reference for merged and joined reads.
    for (genome_len, no_clusters, genome_dir, ref) in index_list:
        if genome_len != 0:
            index_cmd = 'STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s'%(args.threads,genome_dir)
            fasta_file = [file for file in os.listdir(genome_dir) if file.endswith('.fa')][0]
            index_cmd += ' --genomeFastaFiles %s'%os.path.join(genome_dir,fasta_file)
            genomeSAindexNbases = min(14, math.log(genome_len,2)/2 - 1)
            index_cmd += ' --genomeSAindexNbases %i'%genomeSAindexNbases
            print(genome_len)
            genomeChrBinNbits =     min(18, math.log(genome_len/no_clusters, 2))
            index_cmd += ' --genomeChrBinNbits %i' % genomeChrBinNbits
            log = 'making STAR index of %s'%(ref)
            if 'Genome' not in os.listdir(genome_dir):
                run_subprocess([index_cmd], args, log)
    return args

def map_STAR(args):
    """map reads with STAR"""
    for strand in ['watson', 'crick']:
        if strand == 'watson':
            n = 1
        else:
            n = 3
        STAR_index_dir = os.path.join(args.output_dir,'STAR_%s'%(strand))
        cmd = "STAR --runThreadN %s --genomeDir %s"%(args.threads, STAR_index_dir)
        cmd += " --readFilesIn %s " %   vars(args)['%s_r1' % (strand)]
        cmd += " %s" %                  vars(args)['%s_r2' % (strand)]
        cmd += " --outSAMattributes NM MD AS --outSAMtype SAM"
        cmd += " --outFileNamePrefix %s" % (os.path.join(args.output_dir,'%s'%(strand)))
        cmd += " --outReadsUnmapped Fastx" #output of unmapped reads for inspection
        cmd += " --scoreGapATAC -2 --scoreGapNoncan -2"
            #outFilterScoreMinOverLread : float: sam as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads)
            #outFilterMatchNminOverLread: float: same as outFilterScoreMin, but normalized to read length (sum of mates’ lengths for paired-end reads)

            # –outFilterMultimapNmax 1 int: maximum number of loci the read is allowed to map to. Alignments (all of
            # them) will be output only if the read maps to no more loci than this value.
        cmd += " --outFilterMismatchNoverLmax 0.95"
            # TODO: implement --alignEndsType endtoend mapping after joined reads are merged
        cmd += "--outFilterMatchNminOverLread 0.9 --scoreGap -4 " \
                " --alignEndsType Local" \
                " --alignSoftClipAtReferenceEnds No" \
                " --outSAMorder PairedKeepInputOrder" \
                " --outFilterMultimapNmax 1" \
                " --scoreInsOpen -1" \
            #make sure we have a bam file sorted by name
        if args.extraflags:
            cmd += ' %s' % args.extraflags
        log = "run STAR for % strand"%(strand)
        run_subprocess([cmd],args, log)
        log = "write final log of STAR to normal log"
        cmd = "cat %s " % os.path.join(args.output_dir, '%s' % (strand) + 'Log.final.out')
        run_subprocess([cmd], args, log)
    return args


def parse_sam(in_file, out_file, strand):
    """parse sam file and write correct output"""
    out_handle = open(out_file , 'a')
    if strand == 'watson':
        nt = ['C']
    else:
        nt = ['G']
    count = 0
    # print 'Warning, only works for forward mapped reads'
    mismatch = 0
    clip_count_total = 0
    for line in open(in_file, 'r'):
        modulo_line_no = count % 2
        #alternates between 0 and 1
        if line.startswith('@'):
            continue
        split_line = line.rstrip('\n').split('\t')
        #skip read pairs with improper flags.
        #TODO: do this filtering in mark_PCR_duplicates or elsewhere with access to pysam.
        if split_line[1] not in ['0', '99', '147']:
            mismatch += 1
            count += 1
            # continue
        char_count = ''
        clip_count = 0
        for char in split_line[5]:
            if not char.isalpha():
                char_count += char
            elif char == 'S':
                clip_count += int(char_count)
            else:
                char_count = ''
        if clip_count > 6:
            clip_count_total += 1
            count += 1
            # continue
        header = split_line[0].split('|')
        #meth_post list can be present for both R1 and R2 the last Samtools tag added should be the RN:Z: tag, look
        #to the right of this tag only
        meth_pos_list = split_line[0][split_line[0].rindex(':Z:'):].split('|')[1:]
        out_line = [header[0]]
        out_line += split_line[1:9]
        seq = list(split_line[9])
        try:
            meth_pos = [int(n) for n in meth_pos_list[-modulo_line_no].split(',')]
            for n in meth_pos:
                if n >= len(seq):
                    break
                if seq[n] not in ['T','A']:
                    break
                seq[n] = nt[-modulo_line_no]
        except ValueError:
            pass
        out_line += [''.join(seq)]
        out_line += split_line[10:]
        for item in header[1:]:
            if ':' in item and item not in out_line:
                out_line.append(item)
        # out_line += header[3:6]
        out_handle.write('\t'.join(out_line) + '\n')
        count += 1
    print('%s mismatches out of %s' % (mismatch, count))
    print('%s reads out of  %s soft clipped more than 5' % (clip_count_total, count))

def addRG(in_files,args):
    """make header for output bamfile and split in watson and crick"""
    #define readgroup header lines by combining the following

    """
    -
    read group
    ID*
    Unique read group identifier. The value of the ID field is used in the RG tags of alignment records.
    SM*
    Sample (use pool name where a pool is being sequenced)
    LB
    Library
    DS
    Description
    PU
    Platform unit (e.g. lane for Illumina or slide for SOLiD); should be a full, unambiguous identifier
    PI
    Predicted median insert size (maybe different from the actual median insert size)
    CN
    Name of sequencing center producing the read.
    DT
    Date the run was produced (ISO 8601 date or date/time).
    PL
    Platform/technology used to produce the read."""

    with open(args.barcodes,'r') as barcodes:
        sam_out= open(in_files['header'],'a')
        header = barcodes.readline().split('\t')
        for line in barcodes:
            RG = ['@RG']
            split_line = line.split('\t')
            if args.species and 'Species' in header:
                if split_line[(header.index('Species'))] != args.species:
                    continue
            fc = split_line[(header.index('Flowcell'))]
            lane = split_line[(header.index('Lane'))]
            sample = split_line[(header.index('Sample'))]
            RG.append('ID:%s_%s_%s'%(fc,lane,sample))
            RG.append('SM:%s'%(sample))
            RG.append('LB:%s_%s'%(fc,sample))
            RG.append('PL:ILLUMINA\n')
            sam_out.write('\t'.join(RG))
    sam_out.close()
    return in_files


def make_header(args):
    """Make header for watson and crick bam file"""
    header = os.path.join(args.output_dir,'header.sam')
    args.header = header
    header_handle = open(header,'w')
    header_handle.write('@HD\tVN:1.4\n')
    file_sam = open(os.path.join(args.output_dir,'watsonAligned.out.sam'))
    print(file_sam)
    for line in file_sam:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                header_handle.write(line)
            elif not line.startswith('@HD'):
                header_handle.write(line)
        else:
            break
    header_handle.close()
    in_files = {'header':os.path.join(args.output_dir,'header.sam')}
    addRG(in_files, args)
    return args


def bam_output(args):
    """Generate watson and crick output bam file"""

    for strand in ['watson', 'crick']:
        read_sam = os.path.join(args.output_dir,'%sAligned.out.sam' % strand)
        out_sam = tempfile.NamedTemporaryFile(prefix=strand, suffix='.sam', dir=args.output_dir)
        print(out_sam)
        #rewrite sam file merged and joined for watson and crick
        parse_sam(read_sam, out_sam.name, strand)
        #convert to sorted and indexed bam
        cmd = 'cat %s %s |samtools view -@ 4 -Shb |sambamba sort -m 4GB --tmpdir %s -t %s -o %s  /dev/stdin'%(args.header,
                                                                            out_sam.name,args.tmpdir, args.threads,
                                                            os.path.join(args.output_dir,'%s.bam' % strand) )
        log = "make sorted bam file"
        run_subprocess([cmd], args, log)
        out_sam.close()
    return args


def clean(args):
    """delete non-used intermediate files"""
    log =  'removing tmp dir %s ' % (args.tmpdir)
    if args.tmpdir.endswith('STAR'):
        cmd = ['rm -rf %s' % (args.tmpdir)]
        run_subprocess(cmd,args,log)
    log = "remove tmp files from output dir"
    cmd = ['rm -rf %s/*out.sam' % args.output_dir]
    run_subprocess(cmd, args, log)
    cmd = ['rm -rf %s/*out.mate*' % args.output_dir]
    run_subprocess(cmd, args, log)
    cmd = ['rm -rf %s/*Log.out' % args.output_dir]
    run_subprocess(cmd, args, log)
    cmd = ['rm -rf %s/*Log.progress.out' % args.output_dir]
    run_subprocess(cmd, args, log)


def main():
    """main function loop"""
    #1 get command line arguments
    args = parse_args()
    # version = get_version()alignSoftClipAtReferenceEnds
    if __name__ == "__main__":
        log = open(args.log,'w')
    else:
        log = open(args.log, 'a')
    log.write("started run\n")
    #2 make reference genome fo STAR in appropriate directory
    args = index_STAR(args)
    #3 rewrite fastq files to contain
    args = process_reads_watson(args)
    args = process_reads_crick(args)
    #4 map processed reads
    args = map_STAR(args)
    args = make_header(args)
    args = bam_output(args)
    #args = remove_PCR_duplicates(args)
    clean(args)

if __name__ == '__main__':
    main()
