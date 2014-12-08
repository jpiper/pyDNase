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

import argparse, pyDNase
from clint.textui import puts, progress
parser = argparse.ArgumentParser(description='Writes a JSON file of DNase I cuts for regions from a BED file')
parser.add_argument("-w", "--window_size", help="Resize all regions to a specific length",default = 0, type=int)
parser.add_argument("-i",action="store_true", help="Ignores strand information in BED file",default=False)
parser.add_argument("regions", help="BED file of the regions")
parser.add_argument("reads", help="BAM file containing the read data")
parser.add_argument("output", help="filename to write the JSON output to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads)
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
