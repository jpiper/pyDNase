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

import argparse
from clint.textui import progress, puts
import pyDNase

parser = argparse.ArgumentParser(description='Writes two WIG files with the cut information based on the regions in reads BED file and the reads in reads BAM file')
parser.add_argument("-r", "--real",action="store_true", help="Report cuts on the negative strand as positive numbers instead of negative (default: False)",default=False)
parser.add_argument("regions", help="BED file of the regions you want to write wig tracks for")
parser.add_argument("reads", help="The BAM file containing the read data")
parser.add_argument("fw_output", help="Path to write the forward reads wig track to")
parser.add_argument("rev_output", help="Path to write the reverse reads wig track to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads,caching=True)
regions = pyDNase.GenomicIntervalSet(args.regions)
fwigout = open(args.fw_output,"w")
bwigout = open(args.rev_output,"w")

#Prints all the wig values but sorts by chromosome/genomic location first
#TODO: port this most awesome (and hacky) code iteration code to the main API, possibly using a generator expression?
puts("Writing wig tracks...")
for each in progress.bar([item for sublist in sorted(regions.intervals.values()) for item in sorted(sublist, key=lambda peak: peak.startbp)]):
    cuts = reads[each]
    f,r = cuts["+"], cuts["-"]
    #Note that we add 1 to the startbp as WIG is 1-based and the internal logic is 0 based
    print >> fwigout, "fixedStep\tchrom=" + str(each.chromosome) + "\t start="+ str(each.startbp+1) +"\tstep=1"
    for i in f:
        print >> fwigout, i
    print >> bwigout, "fixedStep\tchrom=" + str(each.chromosome) + "\t start="+ str(each.startbp+1) +"\tstep=1"
    for i in r:
        if args.real:
            print >> bwigout, i
        else:
            print >> bwigout, -i

fwigout.close()
bwigout.close()