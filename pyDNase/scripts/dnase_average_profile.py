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
#Required for headless operation
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

parser = argparse.ArgumentParser(description='Plots average profile of DNase activity surrounding a list of regions in a BED file')
parser.add_argument("-w", "--window_size", help="Size of flanking area around centre of the regions to plot (default: 200)",default=200,type=int)
parser.add_argument("-i",action="store_true", help="Ignores any strand information in BED file and plots data relative to reference strand",default=False)
parser.add_argument("regions", help="BED file of the regions you want to generate the average profile for")
parser.add_argument("reads", help="The BAM file containing the DNase-seq data")
parser.add_argument("output", help="filename to write the output to")
args  = parser.parse_args()

xsize   = args.window_size
reads   = pyDNase.BAMHandler(args.reads)
regions = pyDNase.GenomicIntervalSet(args.regions)

#Set all strands to positive if "ignore strands" is enabled
if args.i:
    for each in regions:
        each.strand = "+"

for each in regions:
    each.score = reads.FOS(each)
regions.resizeRegions(xsize)

fw = []
rv = []
for each in regions:
    if reads[each]["+"].sum() and reads[each]["-"].sum():
        fw.append(reads[each]["+"])
        rv.append(reads[each]["-"])
plt.plot(np.mean(fw,axis=0),c="red")
plt.plot(np.mean(rv,axis=0),c="green")

#Pad the axis out reads bit
rcParams['xtick.major.pad'] = 20 
rcParams['ytick.major.pad'] = 20 

#Sort out the X axis
ticks = [0,xsize,xsize*2]
labels = [-xsize,0,xsize]
plt.xticks(ticks, labels)

#Makes ticks only appear on the left hand side
plt.gca().yaxis.set_ticks_position('left')

#Remove top, right, and bottom borders
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
#plt.gca().spines['bottom'].set_visible(False)

plt.gca().tick_params(axis='both', which='major', labelsize=32)

plt.gca().set_ylabel('Average DNase\n Activity',size="36", multialignment='center')
plt.savefig(args.output,bbox_inches='tight')
