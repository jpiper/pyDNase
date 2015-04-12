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

import multiprocessing as mp
import argparse, os
from clint.textui import progress, puts_err
import pyDNase
from pyDNase import footprinting

__version__ = "0.2.0"

parser = argparse.ArgumentParser(description='Footprint the DHSs in a DNase-seq experiment using the Wellington Algorithm.')
parser.add_argument("-b","--bonferroni",action="store_true", help="Performs a bonferroni correction (default: False)",default=False)
parser.add_argument("-sh", "--shoulder-sizes", help="Range of shoulder sizes to try in format \"from,to,step\" (default: 35,36,1)",default="35,36,1",type=str)
parser.add_argument("-fp", "--footprint-sizes", help="Range of footprint sizes to try in format \"from,to,step\" (default: 11,26,2)",default="11,26,2",type=str)
parser.add_argument("-d", "--one_dimension",action="store_true", help="Use Wellington 1D instead of Wellington (default: False)",default=False)
parser.add_argument("-fdr","--FDR_cutoff", help="Write footprints using the FDR selection method at a specific FDR (default: 0.01)",default=0.01,type=float)
parser.add_argument("-fdriter", "--FDR_iterations", help="How many randomisations to use when performing FDR calculations (default: 100)",default=100,type=int)
parser.add_argument("-fdrlimit", "--FDR_limit", help="Minimum p-value to be considered significant for FDR calculation (default: -20)",default=-20,type=int)
parser.add_argument("-pv","--pv_cutoffs", help="Select footprints using a range of pvalue cutoffs (default: -10,-20,-30,-40,-50,-75,-100,-300,-500,-700",default="-10,-20,-30,-40,-50,-75,-100,-300,-500,-700",type=str)
parser.add_argument("-dm","--dont-merge-footprints",action="store_true", help="Disables merging of overlapping footprints (Default: False)",default=False)
parser.add_argument("-o","--output_prefix", help="The prefix for results files (default: <reads.regions>)",default="",type=str)
parser.add_argument("-p", help="Number of processes to use (default: uses all CPUs)",default=0,type=int)
parser.add_argument("regions", help="BED file of the regions you want to footprint")
parser.add_argument("reads", help="The BAM file containing the DNase-seq reads")
parser.add_argument("outputdir", help="A writeable directory to write the results to")
clargs = parser.parse_args()

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
    clargs.shoulder_sizes = xrange_from_string(clargs.shoulder_sizes)
    clargs.footprint_sizes = xrange_from_string(clargs.footprint_sizes)
except ValueError:
    raise RuntimeError("shoulder and footprint sizes must be supplied as from,to,step")

try:
    clargs.pv_cutoffs = map(int,clargs.pv_cutoffs.split(","))
except:
    raise RuntimeError("p-value cutoffs must be supplied as a string of numbers separated by commas")

assert 0 < clargs.FDR_cutoff < 1, "FDR must be between 0 and 1"
assert clargs.FDR_limit < 0, "FDR limit must be less than 0"
assert len([f for f in os.listdir(clargs.outputdir) if f[0] != "."]) == 0, "output directory {0} is not empty!".format(clargs.outputdir)

if not clargs.output_prefix:
    clargs.output_prefix = str(os.path.basename(clargs.reads)) + "." + str(os.path.basename(clargs.regions))

#Load reads and regions
regions = pyDNase.GenomicIntervalSet(clargs.regions)
reads = pyDNase.BAMHandler(clargs.reads,caching=False)

#Create a directory for p-values and WIG output. This /should/ be OS independent
os.makedirs(os.path.join(clargs.outputdir,"p value cutoffs"))
wigout = open(os.path.relpath(clargs.outputdir) + "/" + clargs.output_prefix + ".WellingtonFootprints.wig","w")
fdrout = open(os.path.relpath(clargs.outputdir) + "/" + clargs.output_prefix + ".WellingtonFootprints.FDR.{0}.bed".format(clargs.FDR_cutoff),"w")

#Iterate in chromosome, basepair order
orderedbychr = [item for sublist in sorted(regions.intervals.values()) for item in sorted(sublist, key=lambda peak: peak.startbp)]
puts_err("Calculating footprints...")

if clargs.p:
    CPUs = clargs.p
else:
    CPUs = mp.cpu_count()
max_regions_cached_in_memory = 10 * CPUs
p = mp.Pool(CPUs)

#TODO: How about we store a dictionary of open file handles - or would this cause problems with threading?
def writetodisk(fp):
    #Raw WIG scores
    print >> wigout, "fixedStep\tchrom=" + str(fp.interval.chromosome) + "\t start="+ str(fp.interval.startbp) +"\tstep=1"
    print >> wigout , '\n'.join(map(str, fp.scores))
    #FDR cutoff footprints
    fdr = fp.FDR_value
    if fdr < clargs.FDR_limit:
         for footprint in fp.footprints(withCutoff=fdr,merge=clargs.dont_merge_footprints):
             print >> fdrout, footprint
    #p-value cutoff footprints
    for fpscore in clargs.pv_cutoffs:
        ofile = open(os.path.relpath(os.path.join(clargs.outputdir,"p value cutoffs")) + "/" + clargs.output_prefix + ".WellingtonFootprints.{0}.bed".format(fpscore),"a")
        for footprint in fp.footprints(withCutoff=fpscore):
            print >> ofile, footprint
        ofile.close()

def multiWellington(regions,reads,**kwargs):
    p = mp.Pool(CPUs)
    for i in progress.bar(regions):
        if clargs.one_dimension:
            fp = footprinting.wellington1D(i,reads,**kwargs)
        else:
            fp = footprinting.wellington(i,reads,**kwargs)
        p.apply_async(fp,callback = writetodisk)
        #Hold here while the queue is bigger than the number of reads we're happy to store in memory
        while p._taskqueue.qsize() > max_regions_cached_in_memory:
            pass
    p.close()
    puts_err("Waiting for the last {0} jobs to finish...".format(max_regions_cached_in_memory))
    p.join()

#TODO: Use **args or something similar to pass arguments?
#TODO: Pass FDR Iterations, FDR LimiFDR_cutoff=0.01,FDR_iterations=100
multiWellington(orderedbychr,reads, shoulder_sizes = clargs.shoulder_sizes ,footprint_sizes = clargs.footprint_sizes, FDR_cutoff=clargs.FDR_cutoff,FDR_iterations=clargs.FDR_iterations,bonferroni = clargs.bonferroni)

wigout.close()
