#!/usr/bin/env python
import argparse, pyDNase
from clint.textui import puts, progress
parser = argparse.ArgumentParser(description='Writes a JSON file of DNase I cuts for regions from a BED file')
parser.add_argument("-w", "--window_size", help="Resize all regions to a specific length",default = 0, type=int)
parser.add_argument("-i",action="store_true", help="Ignores strand information in BED file",default=False)
parser.add_argument("-A",action="store_true", help="ATAC-seq mode (default: False)",default=False)
parser.add_argument("regions", help="BED file of the regions")
parser.add_argument("reads", help="BAM file containing the read data")
parser.add_argument("output", help="filename to write the JSON output to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads, ATAC=args.A)
regions = pyDNase.GenomicIntervalSet(args.regions)

if args.i:
    for each in regions:
        each.strand = "+"

if args.window_size:
    puts("Resizing Regions to {0}".format(args.window_size))
    regions.resizeRegions(toSize=args.window_size)

#TODO: this will load everything everything into memory, it's probably worth making the option to write directly to disk

outfile = open(args.output,"w")
outarr = []
puts("Generating JSON output...")
for i in progress.bar(sorted(regions, key = lambda x : x.importorder)):
    cuts = reads[i]
    outarr.append({"location":i.chromosome + ":" + str(i.startbp) + ":" + str(i.endbp), "positive strand cuts":cuts["+"] , "negative strand cuts":cuts["-"]})
puts("Writing JSON to disk...")
outfile.write(str(outarr))
outfile.close()
