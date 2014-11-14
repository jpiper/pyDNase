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
import pyDNase
import numpy as np
import matplotlib as mpl
from clint.textui import progress, puts
#Required for headless operation
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

parser = argparse.ArgumentParser(description='Plots average profile of DNase activity surrounding a list of regions in a BED file')
parser.add_argument("-w", "--window_size", help="Size of flanking area around centre of the regions to plot (default: 100)",default=100,type=int)
parser.add_argument("-bf", "--bias-file", help="Location of the sorted, index",default = None,type=str)
parser.add_argument("-i",action="store_true", help="Ignores any strand information in BED file and plots data relative to reference strand",default=False)
parser.add_argument("-c",action="store_true", help="Combine the strand information into one graph",default=False)
parser.add_argument("-n",action="store_true", help="Normalise cut counts to a fraction peaks",default=False)
parser.add_argument("-b",action="store_true", help="Normalise for cutting bias",default=False)
parser.add_argument("regions", help="BED file of the regions you want to generate the average profile for")
parser.add_argument("reads", help="The BAM file containing the DNase-seq data")
parser.add_argument("output", help="filename to write the output to")
args  = parser.parse_args()

reads   = pyDNase.BAMHandler(args.reads)
if args.b:
    if args.bias_file != None:
        freads   = pyDNase.BAMHandlerWithBias(pyDNase.FASTAHandler(args.bias_file),args.reads)
    else:
        raise ValueError("No FASTA file provided for bias correction!")
regions = pyDNase.GenomicIntervalSet(args.regions)



#Set all strands to positive if "ignore strands" is enabled
if args.i:
    for each in regions:
        each.strand = "+"

puts("Resizing Regions to {0}".format(args.window_size))
regions.resizeRegions(args.window_size)

fw = []
rv = []
puts("Reading Data from BAM file...")
for each in progress.bar(regions):
    if sum(reads[each]["+"]) and sum(reads[each]["-"]):
        if args.b:
            try:
                fw.append(np.divide(reads[each]["+"],freads[each]["+"]))
                rv.append(np.divide(reads[each]["-"],freads[each]["-"]))
            except Exception:
                pass
        else:
            fw.append(reads[each]["+"])
            rv.append(reads[each]["-"])

if args.n:
    fw = [map(float,i)for i in fw]
    rv = [map(float,i) for i in rv]
    fw = [np.divide(np.subtract(i, min(i)), np.subtract(max(i) , min(i))) for i in fw]
    rv = [np.divide(np.subtract(i, min(i)), np.subtract(max(i) , min(i))) for i in rv]

if args.c:
    plt.plot(np.add(np.mean(fw,axis=0),np.mean(rv,axis=0)),c="red")
else:
    plt.plot(np.mean(fw,axis=0),c="red")
    plt.plot(np.mean(rv,axis=0),c="blue")

#Pad the axis out reads bit
rcParams['xtick.major.pad'] = 20 
rcParams['ytick.major.pad'] = 20

#Sort out the X axis ticks
ticks = [0,args.window_size,args.window_size*2]
labels = [-args.window_size,0,args.window_size]
plt.xticks(ticks, labels)

#Make the yaxis start from 0
plt.gca().set_ylim(0)

#Makes ticks only appear on the left hand side
plt.gca().yaxis.set_ticks_position('left')

#Remove top and right borders
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.gca().tick_params(axis='both', which='major', labelsize=28, pad=12)

if args.bias_file:
    plt.gca().set_ylabel('Average DNase Activity\n (Observed/Expected)',size="32", multialignment='center')
else:
    plt.gca().set_ylabel('Average DNase Activity',size="26", multialignment='center')
plt.savefig(args.output,bbox_inches='tight')
