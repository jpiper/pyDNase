#!/usr/bin/env python

# Copyright (C) 2013 Jason Piper - j.piper@warwick.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse, pyDNase, csv, random
import numpy as np
from clint.textui import puts, progress
import warnings

parser = argparse.ArgumentParser(description='Writes a JavaTreeView file based on the regions in reads BED file and the reads in reads BAM file')
parser.add_argument("-w", "--window_size", help="Size of flanking area around centre of the regions to plot (default: 100)",default=100,type=int)
parser.add_argument("-i",action="store_true", help="Ignores strand information in BED file",default=False)
parser.add_argument("-o",action="store_true", help="Orders output the same as the input (default: orders by FOS)",default=False)
parser.add_argument("-a",action="store_true", help="Write absolute cut counts instead strand imbalanced counts",default=False)
parser.add_argument("-n",action="store_true", help="Normalise the cut data for each region between 0 and 1",default=False)
parser.add_argument("-c",action="store_true", help="Disable memory caching (low RAM mode)",default=False)
parser.add_argument("-b",action="store_true", help="Normalise for cutting bias",default=False)
parser.add_argument("-bf", "--bias-file", help="Location of the sorted, index",default = None,type=str)
parser.add_argument("-r",action="store_true", help="Randomise the ordering of the output",default=False)
parser.add_argument("regions", help="BED file of the regions you want to generate the heatmap for")
parser.add_argument("reads", help="The BAM file containing the read data")
parser.add_argument("output", help="filename to write the CSV output to")
args = parser.parse_args()

reads   = pyDNase.BAMHandler(args.reads,caching= not args.c)
if args.b:
    if args.bias_file != None:
        freads   = pyDNase.BAMHandlerWithBias(pyDNase.FASTAHandler(args.bias_file),args.reads,caching= not args.c)
    else:
        raise ValueError("No FASTA file provided for bias correction!")
regions = pyDNase.GenomicIntervalSet(args.regions)

if args.i:
    for each in regions:
        each.strand = "+"

if not args.c:
    puts("Caching reads...")
    for i in progress.bar(regions):
       reads[i]

if args.o:
    sorter = lambda x : x.importorder
else:
    if args.c:
        warnings.warn("You've disabled memory caching, which can cause a ~2x slowdown when sorting by FOS")
    sorter = lambda x : x.score
    puts("Ordering intervals by FOS")
    for i in progress.bar(regions):
        i.score = reads.FOS(i)

puts("Resizing Regions to {0}".format(args.window_size))
regions.resizeRegions(toSize=args.window_size)

outfile = csv.writer(open(args.output,"w"),delimiter="\t")
#Writes the header file
outfile.writerow(["GID", "ID", "NAME"] + [i+1 for i in range(len(reads[regions[0]]["+"]))])

def normalise_cuts(cuts):
    normalised_cuts = np.array(cuts,dtype="float")
    normalised_cuts = ((normalised_cuts - normalised_cuts.min()) / (normalised_cuts.max() - normalised_cuts.min()))
    return normalised_cuts

puts("Writing to JTV...")

if args.r:
    # This is done instead of using a random key to sorted() as it is O(n), whereas sorted() is O(n log n)
    region_iterator = [i for i in regions]
    random.shuffle(region_iterator)
else:
    region_iterator = sorted(regions, key = sorter)

for i in progress.bar(region_iterator):
    cuts = reads[i]
    if args.a:
        newarray = np.add(cuts["+"], cuts["-"])
        if args.n:
            newarray = normalise_cuts(newarray)
    else:
        if args.n:
            cuts["+"] = normalise_cuts(cuts["+"])
            cuts["-"] = normalise_cuts(cuts["-"])

        newarray = np.subtract(cuts["+"], cuts["-"])

    outfile.writerow(["NULL","NULL",i.chromosome + ":" + str(i.startbp) + ":" + str(i.endbp)] + newarray.tolist())