import argparse
import csv
import subprocess
import os

####Parse args
parser = argparse.ArgumentParser(description='Process input files')
parser.add_argument("--r1_in", metavar="reads1", action="store",
                    dest="reads1", help="left-hand fastq file")
parser.add_argument("--r2_in", metavar="reads2", action="store",
                    dest="reads2",
                    help="right-hand fastq file")
parser.add_argument("-b", "--barcodes", metavar="input", action="store",
                    dest="barcode", default="barcodes.tsv",
                    help="input tab separated barcode file")
parser.add_argument("--output-dir", metavar="outputdir", action="store",
                    dest="outputdir", default="",
                    help="Specify output directory, only for galaxy")
parser.add_argument('--tmpdir',
                    help='temporary directory')
args = parser.parse_args()
if args.outputdir:
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
if args.tmpdir:
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    if not os.path.exists(os.path.join(args.outputdir,"clone")):
        os.mkdir(os.path.join(args.outputdir, "clone"))
    if not os.path.exists(os.path.join(args.outputdir,"clone-stacks")):
        os.mkdir(os.path.join(args.outputdir,"clone-stacks"))
if os.path.exists(os.path.join(args.outputdir, "demultiplex.log")):
    os.remove(os.path.join(args.outputdir, "demultiplex.log"))

def run_subprocess(cmd, args, log_message):
    "Run subprocess under standardized settings"
    # force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(os.path.join(args.outputdir, "demultiplex.log"), 'a') as log:
        log.write("now starting:\t%s\n" % log_message)
        log.write('running:\t%s\n' % (' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode().replace('\r', '\n')
        stderr = stderr.decode().replace('\r', '\n')
        if stdout:
            log.write('stdout:\n%s\n' % stdout)
        if stderr:
            log.write('stderr:\n%s\n' % stderr)
        log.write('finished:\t%s\n\n' % log_message)
    return 0



# Get info from barcode file
bc_file = open(args.barcode, 'r')
header = bc_file.readline().rstrip("\n")
bc_dict = {}
line2 = (bc_file.readlines()[1])
head = header.split(sep="\t")
line2split = line2.split(sep="\t")
nCol_bc = range(len(head))

for i in nCol_bc:
    bc_dict[head[i]] = [i]

# Make barcode as stacks input #################
with open(os.path.join(args.outputdir, "barcode_stacks.tsv"), 'w', newline='') as out_file:
    ##TODO make proper tmp file
    tsv_writer = csv.writer(out_file, delimiter='\t')
    for strand in ["watson", "crick"]:
        bc_file = open(args.barcode, 'r')
        lines = bc_file.readlines()[1:]
        for line in lines:
            if strand == 'watson':
                lineout = line.split()
                valuesout = []
                valuesout.append(lineout[int(bc_dict["Barcode_R1"][0])] + "T")
                valuesout.append(lineout[int(bc_dict["Barcode_R2"][0])] + "C")
                valuesout.append(lineout[int(bc_dict["Sample"][0])] + "-Watson")
                tsv_writer.writerow(valuesout)
            else:
                lineout = line.split()
                valuesout = []
                valuesout.append(lineout[bc_dict["Barcode_R1"][0]] + "C")
                valuesout.append(lineout[bc_dict["Barcode_R2"][0]] + "T")
                valuesout.append(lineout[bc_dict["Sample"][0]] + "-Crick")
                tsv_writer.writerow(valuesout)
out_file.close()

# Run stacks ###############
# Run clone filter

if int(line2split[bc_dict["Wobble_R1"][0]]) is not 0:
    cmd = "clone_filter -1 %s -2 %s " % (args.reads1, args.reads2)
    cmd += "-o %s --inline_inline " % os.path.join(args.outputdir, "clone")
    cmd += "-igzfastq --oligo_len_1 %s --oligo_len_2 %s " % (line2split[bc_dict["Wobble_R1"][0]],
                                                            line2split[bc_dict["Wobble_R2"][0]])
    log = "Run clone_filter from stacks"
    print(log)
    run_subprocess([cmd], args, log)
# Run remove radtags

if int(line2split[bc_dict["Wobble_R1"][0]]) is 0:
    read1 = args.reads1
    read2 = args.reads2
else:
    cloneList = os.listdir(os.path.join(args.outputdir, "clone"))  # find files in the clone directory
    for clonefile in cloneList:  # "Grep" non-discarded
        if ".1.fq.gz" in clonefile:
            read1 = os.path.join(args.outputdir, "clone",clonefile)
        if ".2.fq.gz" in clonefile:
            read2 = os.path.join(args.outputdir, "clone", clonefile)


cmd = "process_radtags -1 %s -2 %s " % (read1,
                                        read2)
cmd += "-b %s -o %s " % (os.path.join(args.outputdir, "barcode_stacks.tsv"), os.path.join(args.outputdir, "clone-stacks"))
cmd += "-r -D --inline_inline --barcode_dist_2 0 "
cmd += "--renz_1 %s --renz_2 %s --retain_header --barcode_dist_2 0" % (line2split[bc_dict["ENZ_R1"][0]],
                                                     line2split[bc_dict["ENZ_R2"][0]])
log = "Run process_radtags from stacks"
print(log)
run_subprocess([cmd], args, log)


