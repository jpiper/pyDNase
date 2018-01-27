import argparse
import pyDNase
from clint.textui import progress

parser = argparse.ArgumentParser(description='writes a BED file with the FOS for the interval specified as the score')
parser.add_argument("-A",action="store_true", help="ATAC-seq mode (default: False)",default=False)
parser.add_argument("regions", help="BED file of the regions you want to generate the average profile for")
parser.add_argument("reads", help="The BAM file containing the DNase-seq data")
parser.add_argument("output", help="filename to write the output to")
args  = parser.parse_args()

reads   = pyDNase.BAMHandler(args.reads, ATAC=args.A)
regions = pyDNase.GenomicIntervalSet(args.regions)

outfile = open(args.output,"w")
for i in progress.bar(regions):
    i.score = reads.FOS(i)
    print(i, file=outfile)
