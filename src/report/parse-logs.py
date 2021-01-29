import argparse
import os 
import re
def parse_args():
    #"""Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Parse input for MultiQC, e.g. python ./src/report/parse-logs.py -p ./output/output_demultiplex/clone-stacks/process_radtags.clone.log -m ./output/output_denovo/make_reference.log -d ./output/output_demultiplex/demultiplex.log -o output/log')
    parser.add_argument('-p', '--input_stacks_bc', type=str, default='process_radtags.clone.log',
                        help='process_radtags.clone.log')
    parser.add_argument('-d', '--input_stacks', type=str, default='demultiplex.log',
                        help='demultiplex.log')
    parser.add_argument('-m', '--input_denovo', type=str, default='make_reference.log',
                        help='make_reference.log')
    parser.add_argument('-o', '--output_dir', type=str, default='log',
                        help='output directory')
    args = parser.parse_args()
    return args
​
def parse_stacks_bc(args):
    outdir = args.output_dir
    output_stacks_bc = open(os.path.join(outdir, "stacks_denovo_barcodes.log"), "w")
    input_stacks_bc = open(args.input_stacks_bc)
    for line in input_stacks_bc:
        if line.startswith("Sequences not recorded"):
            for _ in range(21):
                # first put next() instead of readline(), next() is python 2.7!
                output_stacks_bc.write(input_stacks_bc.readline())
​
    output_stacks_bc.close()
​
def parse_stacks(args):
    outdir = args.output_dir
    output = open(os.path.join(outdir, "demultiplexing.log"), "w")
    output.write("name, count, %\n")
    output_clones = open(os.path.join(outdir, "clones.log"), "w")
    input_stacks = open(args.input_stacks)
    for line in input_stacks:
        #if "total sequences" in line:
         #   output.write(line.lstrip())
        if "ambiguous barcode drops (" in line:
            string = line.lstrip()
            list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
            output.write("ambiguous barcode drops, " + ", ".join(list) + '\n')
        elif "low quality read drops (" in line:
            string = line.lstrip()
            list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
            output.write("low quality read drops, " + ", ".join(list) + '\n')        
        elif "ambiguous RAD-Tag drops (" in line:
            string = line.lstrip()
            list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
            output.write("ambiguous RAD-Tag drops, " + ", ".join(list) + '\n')        
        elif "retained reads (" in line:
            string = line.lstrip()
            list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
            output.write("retained reads, " + ", ".join(list) + '\n')        
        elif "Calculating the distribution of cloned read pairs..." in line:
            for _ in range(1):
                #output_clones.write(input_stacks.readline())
                string = input_stacks.readline()
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_clones.write("input read pairs, output read pairs, discarded read pairs, % clone reads\n" + ", ".join(list)) 
​
    output.close()   
    output_clones.close()
​
def parse_denovo(args):
    outdir = args.output_dir
    output_denovo = open(os.path.join(outdir, "denovo.log"), "w")
    output_merge = open(os.path.join(outdir,"merge.log"), "w")
    output_merge.write("name, count, total count, %\n")
    output_denovo.write("total clusters, min size, max size, avg. size, singletons, singletons % seq, singletons % clusters\n") 
    count = 0
    count2 = 0
    count3 = 0
    clusters = []
    for line in open(args.input_denovo):
        if "Assembled reads ." in line:
            count += 1
            if count == 1:
                strand = 'crick'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Assembled reads " + strand + ", " + ", ".join(list) + '\n')
            elif count == 2:
                strand = 'watson'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Assembled reads " + strand + ", " + ", ".join(list) + '\n')
        elif "Discarded reads ." in line:
            count2 += 1
            if count2 == 1:
                strand = 'crick'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Discarded reads " + strand + ", " + ", ".join(list) + '\n')
            elif count2 == 2:
                strand = 'watson'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Discarded reads " + strand + ", " + ", ".join(list) + '\n')
        elif "Not assembled reads ." in line:
            count3 += 1
            if count3 == 1:
                strand = 'crick'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Not assembled reads " + strand + ", " + ", ".join(list) + '\n')
            elif count3 == 2:
                strand = 'watson'
                string = line.replace(',', '')
                list = [str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string)
                output_merge.write("Not assembled reads " + strand + ", " + ", ".join(list) + '\n')
        elif "Clusters:" in line:
                string = line.replace(',', '')
                clusters.extend([str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string))
        elif "Singletons:" in line:
                string = line.replace(',', '')
                clusters.extend([str(s) for s in string.split() if s.isdigit()] + re.findall(r'\d+\.\d+', string))
    #print(clusters)
    output_denovo.write(", ".join(clusters))
​
​
    output_denovo.close()
    output_merge.close()
​
​
def main():
    args = parse_args()
    parse_stacks_bc(args)
    parse_stacks(args)
    parse_denovo(args)
​
if __name__ == "__main__":
    main()
