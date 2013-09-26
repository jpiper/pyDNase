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

import argparse, pyDNase, csv
from clint.textui import puts, progress
parser = argparse.ArgumentParser(description='Writes a JavaTreeView file based on the regions in reads BED file and the reads in reads BAM file')
parser.add_argument("-w", "--window_size", help="Size of flanking area around centre of the regions to plot (default: 100)",default=100,type=int)
parser.add_argument("-i",action="store_true", help="Ignores strand information in BED file",default=False)
parser.add_argument("-o",action="store_true", help="Orders output the same as the input (default: orders by FOS)",default=False)
parser.add_argument("-a",action="store_true", help="Write absolute cut counts instead strand imbalanced counts",default=False)
parser.add_argument("regions", help="BED file of the regions you want to generate the heatmap for")
parser.add_argument("reads", help="The BAM file containing the read data")
parser.add_argument("output", help="filename to write the CSV output to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads)
regions = pyDNase.GenomicIntervalSet(args.regions)

if args.i:
    for each in regions:
        each.strand = "+"

puts("Caching reads...")
for i in progress.bar(regions):
   reads[i]

if args.o:
    sorter = lambda x : x.importorder
else:
    sorter = lambda x : x.score
    puts("Ordering intervals by FOS")
    for i in progress.bar(regions):
        i.score = reads.FOS(i)

puts("Resizing Regions to {}".format(args.window_size))
regions.resizeRegions(toSize=args.window_size)

outfile = csv.writer(open(args.output,"w"),delimiter="\t")
#Writes the header file
outfile.writerow(["GID", "ID", "NAME"] + [i+1 for i in range(len(reads[regions[0]]["+"]))])


puts("Writing to JTV...")
for i in progress.bar(sorted(regions, key = sorter)):
    cuts = reads[i]
    if args.a:
        newarray = (cuts["+"] + cuts["-"]).tolist()
    else:
        newarray = (cuts["+"] - cuts["-"]).tolist()
    outfile.writerow(["NULL","NULL",i.chromosome + ":" + str(i.startbp) + ":" + str(i.endbp)] + newarray)