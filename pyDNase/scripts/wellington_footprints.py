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

import argparse, os
from scipy import stats
from clint.textui import progress,puts
import pyDNase
from pyDNase import footprinting

__version__ = "0.1.0"

parser = argparse.ArgumentParser(description='Footprint the DHSs in a DNase-seq experiment using the Wellington Algorithm.')
parser.add_argument("-b","--bonferroni",action="store_true", help="Performs a bonferroni correction (default: False)",default=False)
parser.add_argument("-sh", "--shoulder-sizes", help="Range of shoulder sizes to try in format \"from,to,step\" (default: 35,36,1)",default="35,36,1",type=str)
parser.add_argument("-fp", "--footprint-sizes", help="Range of footprint sizes to try in format \"from,to,step\" (default: 11,26,2)",default="11,26,2",type=str)
parser.add_argument("-d", "--one-dimension",action="store_true", help="Use Wellington 1D instead of Wellington (default: False)",default=False)
parser.add_argument("-fdr","--FDR_cutoff", help="Write footprints using the FDR selection method at a specific FDR (default: 0.01)",default=0.01,type=float)
parser.add_argument("-fdriter", "--FDR-iterations", help="How many randomisations to use when performing FDR calculations (default: 100)",default=100,type=int)
parser.add_argument("-fdrlimit", "--FDR-limit", help="Minimum p-value to be considered significant for FDR calculation (default: -20)",default=-20,type=int)
parser.add_argument("-pv","--pv_cutoffs", help="Select footprints using a range of pvalue cutoffs (default: -10,-20,-30,-40,-50,-75,-100,-300,-500,-700",default="-10,-20,-30,-40,-50,-75,-100,-300,-500,-700",type=str) #map(int,"1,2,3".split(","))
parser.add_argument("-dm","--dont-merge-footprints",action="store_true", help="Disables merging of overlapping footprints (Default: False)",default=False)
parser.add_argument("-o","--output_prefix", help="The prefix for results files (default: <reads.regions>)",default="",type=str)
parser.add_argument("regions", help="BED file of the regions you want to footprint")
parser.add_argument("reads", help="The BAM file containing the DNase-seq reads")
parser.add_argument("outputdir", help="A writeable directory to write the results to")
args = parser.parse_args()

#Sanity check parameters from the user

def xrange_from_string(range_string):
    try:
        range_string = map(int,range_string.split(","))
        range_string = range(range_string[0],range_string[1],range_string[2])
        assert len(range_string) > 0
        return range_string
    except:
        raise ValueError

try:
    args.shoulder_sizes = xrange_from_string(args.shoulder_sizes)
    args.footprint_sizes = xrange_from_string(args.footprint_sizes)
except ValueError:
    raise RuntimeError("shoulder and footprint sizes must be supplied as from,to,step")

try:
    args.pv_cutoffs = map(int,args.pv_cutoffs.split(","))
except:
    raise RuntimeError("p-value cutoffs must be supplied as a string of numbers separated by commas")

assert 0 < args.FDR_cutoff < 1, "FDR must be between 0 and 1"
assert args.FDR_limit < 0, "FDR limit must be less than 0"

#Checks that the directories are empty (ignores hidden files/folders)
assert len([f for f in os.listdir(args.outputdir) if f[0] != "."]) == 0, "output directory {0} is not empty!".format(args.outputdir)

if not args.output_prefix:
    args.output_prefix = str(os.path.basename(args.reads)) + "." + str(os.path.basename(args.regions))

#Load reads and regions
regions = pyDNase.GenomicIntervalSet(args.regions)
reads = pyDNase.BAMHandler(args.reads,caching=False)

#Create a directory for p-values and WIG output. This /should/ be OS independent
os.makedirs(os.path.join(args.outputdir,"p value cutoffs"))
wigout = open(os.path.relpath(args.outputdir) + "/" + args.output_prefix + ".WellingtonFootprints.wig","w")
fdrout = open(os.path.relpath(args.outputdir) + "/" + args.output_prefix + ".WellingtonFootprints.FDR.{0}.bed".format(args.FDR_cutoff),"w")

#Iterate in chromosome, basepair order
orderedbychr = [item for sublist in sorted(regions.intervals.values()) for item in sorted(sublist, key=lambda peak: peak.startbp)]
puts("Calculating footprints...")
for each in progress.bar(orderedbychr):
    #Calculate footprint scores (1D or 2D)
    #TODO: put args here.
    if args.one_dimension:
        fp = footprinting.wellington1D(each, reads, shoulder_sizes = args.shoulder_sizes ,footprint_sizes = args.footprint_sizes, bonferroni = args.bonferroni)
    else:
        fp = footprinting.wellington(each, reads, shoulder_sizes = args.shoulder_sizes ,footprint_sizes = args.footprint_sizes, bonferroni = args.bonferroni)

    #Write fpscores to WIG
    print >> wigout, "fixedStep\tchrom=" + str(fp.interval.chromosome) + "\t start="+ str(fp.interval.startbp) +"\tstep=1"
    for i in fp.scores:
        print >> wigout, i

    #FDR footprints
    fdr = stats.scoreatpercentile([a for a in fp.calculate(reads,FDR=True, shoulder_sizes = args.shoulder_sizes ,footprint_sizes = args.footprint_sizes, bonferroni = args.bonferroni)[0].tolist() for i in range(args.FDR_iterations)],int(args.FDR_cutoff*100))
    if fdr < args.FDR_limit:
        for footprint in fp.footprints(withCutoff=fdr,merge=not args.dont_merge_footprints):
            print >> fdrout, footprint

    #p-value cutoff footprints
    for fpscore in args.pv_cutoffs:
        ofile = open(os.path.relpath(os.path.join(args.outputdir,"p value cutoffs")) + "/" + args.output_prefix + ".WellingtonFootprints.{0}.bed".format(fpscore),"a")
        for footprint in fp.footprints(withCutoff=fpscore):
            print >> ofile, footprint
        ofile.close()
wigout.close()
