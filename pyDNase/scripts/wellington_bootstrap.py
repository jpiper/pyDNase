#!/usr/bin/env python

import pyDNase
from pyDNase import footprinting
import pyDNase.footprinting.WellingtonC
import numpy as np
from clint.textui import progress
import argparse, random

parser = argparse.ArgumentParser(description='Scores Differential Footprints using Wellington-Bootstrap')
parser.add_argument("treatment",help="BAM file for treatment")
parser.add_argument("control",help="BAM file for control")
parser.add_argument("bedsites",help="BED file of genomic locations to search in")
parser.add_argument("treatment_only_output",help="File to write treatment specific fooprints scores to")
parser.add_argument("control_only_output",help="File to write control specific footprint scores to")
args  = parser.parse_args()
import multiprocessing as mp


#Treatment
reads2 = pyDNase.BAMHandler(args.treatment,caching=0)
#Control
reads1  = pyDNase.BAMHandler(args.control,caching=0)
#Regions of Interest
regions = pyDNase.GenomicIntervalSet(args.bedsites)
#output
t_ofile = open(args.treatment_only_output,"w",buffering=1)
c_ofile = open(args.control_only_output,"w",buffering=1)


shuffle_1 = 0
MIN_SCORE = -20
self_scoring = 0


class Diffwell(footprinting.wellington):
    def __init__(self, reads2, *args, **kwargs):
        super(Diffwell, self).__init__(*args, **kwargs)
        self.reads2 = reads2[self.interval]
       # self.bias_data     = BiasCalculator()

    def footprints(self, withCutoff, merge = 0,select = 1):
        """
        This returns reads GenomicIntervalSet with the intervals retrieved below the specific cutoff applied to the selected data
        """
        #This find the positions of all the ranges below the cutoff using reads new method

        ranges = []

        tempMLE, templogProb = np.copy(self.lengths), np.copy(self.scores)

        while templogProb.min() < withCutoff:
            minimapos = templogProb.argmin()
            minimafplen = tempMLE[minimapos]
            minimaphalffplen = int(minimafplen)/2
            lbound = max(minimapos-(minimaphalffplen),0)
            rbound = min(minimapos+(minimaphalffplen),len(templogProb))
            ranges.append((lbound,rbound,templogProb.min(),minimafplen))
            templogProb[max(lbound-self.shoulder_sizes[-1],0):min(rbound+self.shoulder_sizes[-1],len(templogProb))] = 1

        returnSet = []
        if ranges:
            if merge:
                ranges = sorted(ranges)
                merged_ranges = [ranges[0]]
                for c, d, e in ranges[1:]:
                    a, b, f = merged_ranges[-1]
                    if c<=b<d:
                        merged_ranges[-1] = a, d , min(e,f)
                    elif b<c<d:
                        merged_ranges.append((c,d,e))
            elif select:
                merged_ranges = []
                while len(ranges):
                    #Find best score
                    sorted(ranges,key=lambda x:-x[2])
                    #Take the last value
                    best = ranges.pop()
                    merged_ranges.append(best)
                    #Check for overlapping regions and remove
                    new_ranges = []
                    for c, d, e, f in ranges:
                        if not c<=best[1]<=d:
                            new_ranges.append([c, d, e,f])
                    ranges = new_ranges
            else:
                 merged_ranges = ranges
            #Creates reads GenomicIntervalSet and adds the footprints to them
            for i in merged_ranges:
                returnSet.append(((i[0] + i[1])/2, i[3]))
        return returnSet

    def findDiffFP(self):
        cuts = self.reads
        forwardArray, backwardArray     = cuts["+"], cuts["-"]
        if self_scoring:
            forwardArray2, backwardArray2 = forwardArray, backwardArray
        elif shuffle_1:
            forwardArray2, backwardArray2 = cuts["+"][:], cuts["-"][:]
            random.shuffle(forwardArray2)
            random.shuffle(backwardArray2)
        else:
            cuts2 = self.reads2
            forwardArray2, backwardArray2 = cuts2["+"], cuts2["-"]

        # Adjust the FDR threshold to a minimum of MIN_SCORE
        threshold = min(self.FDR_value, MIN_SCORE)

        offsets = self.footprints(threshold,select=1)

        # This is where the magic happens
        best_probabilities,best_footprintsizes= pyDNase.footprinting.WellingtonC.diff_calculate(forwardArray,backwardArray,forwardArray2,backwardArray2,[i[1] for i in offsets],[i[0] for i in offsets],threshold)

        rset = []

        for i in offsets:
            middle = self.interval.startbp + i[0]
            fp_halfsize = (best_footprintsizes[i[0]] // 2)
            left = middle - fp_halfsize
            right = middle + fp_halfsize
            mmscore = best_probabilities[i[0]]
            heyo = pyDNase.GenomicInterval(self.interval.chromosome,left,right,score=mmscore)
            rset.append(heyo)
        return rset

    def __call__(self):
        results = None
        if min(self.scores) < -20: #Previously -20
            if min(self.scores) < self.FDR_value:
                results = self.findDiffFP()
        self.results = results

        return(self)


CPUs = mp.cpu_count()

# This roughly scales at aabout 450mb per 300 regions held in memory
max_regions_cached_in_memory = 50 * CPUs
p = mp.Pool(CPUs)

print "Performing Differential Footprinting..."

def write_treat_to_disk(item):
    if item.results:
        for i in item.results:
            print >> t_ofile, i

def write_control_to_disk(item):
    if item.results:
        for i in item.results:
            print >> c_ofile, i

for i in progress.bar(regions):
    # Make sure the interval is actually big enough to footprint to begin with
    if len(i) < 120:
        i.startbp -= 60
        i.endbp += 60
    # Perform both comparisons - A against B and B against A
    fp  = Diffwell(reads2=reads1, interval=i, reads=reads2)
    fp2 = Diffwell(reads2=reads2, interval=i, reads=reads1)
    p.apply_async(fp, callback = write_treat_to_disk)
    p.apply_async(fp2, callback = write_control_to_disk)
    #Hold here while the queue is bigger than the number of reads we're happy to store in memory
    while p._taskqueue.qsize() > max_regions_cached_in_memory:
        pass
