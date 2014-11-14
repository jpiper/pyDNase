#!/usr/bin/env python

# Copyright (C) 2014 Jason Piper - j.piper@warwick.ac.uk
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
import pyDNase
import math
from clint.textui import progress, puts

parser = argparse.ArgumentParser(description='Annotates a set of DHSs with the dDHS score (He et al. 2012)')
parser.add_argument("-l", action="store_true", help="low RAM mode (disables caching) (Default: False)",default=False)
parser.add_argument("regions", help="The set of BED files you wish to annotate with dDHS scores")
parser.add_argument("treat_dhs", help="The DHSs belonging to the Treatment")
parser.add_argument("control_dhs", help="The DHSs belonging to the control")
parser.add_argument("reads_treat", help="The BAM file containing the Treatment DNase-seq data")
parser.add_argument("reads_control", help="The BAM file containing the Control DNase-seq data")
parser.add_argument("output", help="filename to write the output to")
args  = parser.parse_args()

reads_treat   = pyDNase.BAMHandler(args.reads_treat, caching = not args.l)
reads_control = pyDNase.BAMHandler(args.reads_control, caching = not args.l)
treat_dhs     = pyDNase.GenomicIntervalSet(args.treat_dhs)
control_dhs   = pyDNase.GenomicIntervalSet(args.control_dhs)
regions       = pyDNase.GenomicIntervalSet(args.regions)

treat_total_cuts   = 0
control_total_cuts = 0
treat_base_pairs   = 0
control_base_pairs = 0

puts("Calculating enrichment for Treatment")
for i in progress.bar(treat_dhs):
    treat_total_cuts += sum([sum(j) for j in reads_treat[i].values()])
    treat_base_pairs += len(i)

puts("Calculating enrichment for Control")
for i in progress.bar(control_dhs):
    control_total_cuts += sum([sum(j) for j in reads_control[i].values()])
    control_base_pairs += len(i)

treat_base_pairs   = float(treat_base_pairs)
control_base_pairs = float(control_base_pairs)
treat_total_cuts   = float(treat_total_cuts)
control_total_cuts = float(control_total_cuts)

ofile = open(args.output,"w")

puts("Calculating dDHS scores...")
for i in progress.bar(regions):
    treat_cuts   =  sum([sum(j) for j in reads_treat[i].values()])
    control_cuts =  sum([sum(j) for j in reads_control[i].values()])
    LHS = math.sqrt(  treat_cuts /(  treat_total_cuts / treat_base_pairs))
    RHS = math.sqrt(control_cuts /(control_total_cuts / control_base_pairs))
    i.score = LHS - RHS
    print >> ofile, i
ofile.close()
