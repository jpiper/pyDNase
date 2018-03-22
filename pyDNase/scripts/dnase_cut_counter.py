#!/usr/bin/env python
import argparse
import pyDNase
from clint.textui import progress, puts

parser = argparse.ArgumentParser(description='Annotates a BED file with with the number of DNase cuts in it')
parser.add_argument("-A",action="store_true", help="ATAC-seq mode (default: False)",default=False)
parser.add_argument("regions", help="BED file")
parser.add_argument("reads", help="The BAM file containing the DNase-seq data")
parser.add_argument("output", help="filename to write the output to")
args  = parser.parse_args()

reads  = pyDNase.BAMHandler(args.reads, caching=0, ATAC=args.A)
regions = pyDNase.GenomicIntervalSet(args.regions)
ofile = open(args.output,"w")

puts("Writing output...")
for i in progress.bar(regions):
    i.score   = sum([sum(j) for j in list(reads[i].values())])
    print(i, file=ofile)
