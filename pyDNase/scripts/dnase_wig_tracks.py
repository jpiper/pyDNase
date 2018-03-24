#!/usr/bin/env python
import argparse
from clint.textui import progress, puts
import pyDNase

parser = argparse.ArgumentParser(description='Writes two WIG files with the cut information based on the regions in reads BED file and the reads in reads BAM file')
parser.add_argument("-r", "--real",action="store_true", help="Report cuts on the negative strand as positive numbers instead of negative (default: False)",default=False)
parser.add_argument("-A",action="store_true", help="ATAC-seq mode (default: False)",default=False)
parser.add_argument("regions", help="BED file of the regions you want to write wig tracks for")
parser.add_argument("reads", help="The BAM file containing the read data")
parser.add_argument("fw_output", help="Path to write the forward reads wig track to")
parser.add_argument("rev_output", help="Path to write the reverse reads wig track to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads, caching=True, ATAC=args.A)
regions = pyDNase.GenomicIntervalSet(args.regions)
fwigout = open(args.fw_output,"w")
bwigout = open(args.rev_output,"w")

#Required for UCSC upload
print("track type=wiggle_0", file=fwigout)
print("track type=wiggle_0", file=bwigout)

#Prints all the wig values but sorts by chromosome/genomic location first
#TODO: port this most awesome (and hacky) code iteration code to the main API, possibly using a generator expression?
puts("Writing wig tracks...")
for each in progress.bar([item for sublist in sorted(regions.intervals.values()) for item in sorted(sublist, key=lambda peak: peak.startbp)]):
    cuts = reads[each]
    f,r = cuts["+"], cuts["-"]
    #Note that we add 1 to the startbp as WIG is 1-based and the internal logic is 0 based
    print("fixedStep\tchrom=" + str(each.chromosome) + "\t start="+ str(each.startbp+1) +"\tstep=1", file=fwigout)
    for i in f:
        print(i, file=fwigout)
    print("fixedStep\tchrom=" + str(each.chromosome) + "\t start="+ str(each.startbp+1) +"\tstep=1", file=bwigout)
    for i in r:
        if args.real:
            print(i, file=bwigout)
        else:
            print(-i, file=bwigout)

fwigout.close()
bwigout.close()
