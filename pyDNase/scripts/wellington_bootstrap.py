#!/usr/bin/env python

# Copyright (C) 2015 Jason Piper - j.piper@warwick.ac.uk
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

import pyDNase, pyDNase.footprinting
import numpy as np
from clint.textui import progress
import multiprocessing as mp
import argparse

__version__ = "0.1.0"

def write_treat_to_disk(item):
    if item.results:
        for i in item.results:
            print >> treatment_output, i


def write_control_to_disk(item):
    if item.results:
        for i in item.results:
            print >> control_output, i


def xrange_from_string(range_string):
    try:
        range_string = map(int, range_string.split(","))
        range_string = range(range_string[0], range_string[1], range_string[2])
        assert len(range_string) > 0
        return range_string
    except:
        raise ValueError


class Diffwell(pyDNase.footprinting.wellington):
    def __init__(self, reads2, min_score,  *args, **kwargs):
        super(Diffwell, self).__init__(*args, **kwargs)
        self.MIN_SCORE = min_score
        self.reads2 = reads2[self.interval]

    def footprints(self, withCutoff = -20, merge = 1):
        """
        This returns reads GenomicIntervalSet with the intervals retrieved
        below the specific cutoff applied to the selected data
        """

        ranges = []

        tempMLE, templogProb = np.copy(self.lengths), np.copy(self.scores)

        while templogProb.min() < withCutoff:
            minimapos = templogProb.argmin()
            minimafplen = tempMLE[minimapos]
            minimaphalffplen = minimafplen//2
            lbound = max(minimapos-minimaphalffplen, 0)
            rbound = min(minimapos+minimaphalffplen, len(templogProb))
            ranges.append((lbound, rbound, templogProb.min(), minimafplen))
            templogProb[max(lbound-self.shoulder_sizes[-1], 0):min(rbound+self.shoulder_sizes[-1], len(templogProb))] = 1

        return_set = []
        if ranges:
            merged_ranges = []
            while len(ranges):
                # Find best score
                sorted(ranges, key=lambda x: -x[2])
                # Take the last value
                best = ranges.pop()
                merged_ranges.append(best)
                # Check for overlapping regions and remove
                new_ranges = []
                for c, d, e, f in ranges:
                    if not c <= best[1] <= d:
                        new_ranges.append([c, d, e, f])
                ranges = new_ranges
            else:
                 merged_ranges = ranges
            # Creates reads GenomicIntervalSet and adds the footprints to them
            for i in merged_ranges:
                return_set.append(((i[0] + i[1])/2, i[3]))
        return return_set

    def findDiffFP(self):
        cuts = self.reads
        forwardArray, backwardArray = cuts["+"], cuts["-"]

        cuts2 = self.reads2
        forwardArray2, backwardArray2 = cuts2["+"], cuts2["-"]

        # Adjust the FDR threshold to a minimum of withCutoff
        threshold = min(self.FDR_value, self.MIN_SCORE)

        # Find the footprints at this threshold
        offsets = self.footprints(threshold)

        # Work out the bootstrap scores for these footprints using the other data set
        best_probabilities, best_footprintsizes = pyDNase.footprinting.WellingtonC.diff_calculate(forwardArray,
                                                                                                  backwardArray,
                                                                                                  forwardArray2,
                                                                                                  backwardArray2,
                                                                                                  [i[1] for i in offsets],
                                                                                                  [i[0] for i in offsets],
                                                                                                  threshold)

        result_intervals = []

        for i in offsets:
            middle = self.interval.startbp + i[0]
            fp_halfsize = (best_footprintsizes[i[0]] // 2)
            left = middle - fp_halfsize
            right = middle + fp_halfsize
            ml_score = best_probabilities[i[0]]
            result = pyDNase.GenomicInterval(self.interval.chromosome, left, right, score=ml_score)
            result_intervals.append(result)
        return result_intervals

    def __call__(self):
        results = None
        # this is where the first round of footprinting is actually called, as self.scores invoked the footprinting
        if min(self.scores) < self.MIN_SCORE:
            if min(self.scores) < self.FDR_value:
                results = self.findDiffFP()
        self.results = results
        return self



parser = argparse.ArgumentParser(description='Scores Differential Footprints using Wellington-Bootstrap.')
parser.add_argument("-fp", "--footprint-sizes",
                    help="Range of footprint sizes to try in format \"from,to,step\" (default: 11,26,2)",
                    default="11,26,2",
                    type=str)
parser.add_argument("-fdr","--FDR_cutoff",
                    help="Detect footprints using the FDR selection method at a specific FDR (default: 0.01)",
                    default=0.01,
                    type=float)
parser.add_argument("-fdriter", "--FDR_iterations",
                    help="How many randomisations to use when performing FDR calculations (default: 100)",
                    default=100,
                    type=int)
parser.add_argument("-fdrlimit", "--FDR_limit",
                    help="Minimum p-value to be considered significant for FDR calculation (default: -20)",
                    default=-20,
                    type=int)
parser.add_argument("-p", "--processes", help="Number of processes to use (default: uses all CPUs)",
                    default=0,
                    type=int)
parser.add_argument("treatment_bam", help="BAM file for treatment")
parser.add_argument("control_bam", help="BAM file for control")
parser.add_argument("bedsites", help="BED file of genomic locations to search in")
parser.add_argument("treatment_only_output", help="File to write treatment specific fooprints scores to")
parser.add_argument("control_only_output", help="File to write control specific footprint scores to")
args = parser.parse_args()

# Sanity check parameters from the user

try:
    args.footprint_sizes = xrange_from_string(args.footprint_sizes)
except ValueError:
    raise RuntimeError("Footprint sizes must be supplied as from,to,step")

assert 0 < args.FDR_cutoff < 1, "FDR must be between 0 and 1"
assert args.FDR_limit < 0, "FDR limit must be less than 0"

# Treatment
reads2 = pyDNase.BAMHandler(args.treatment_bam, caching=0)
# Control
reads1 = pyDNase.BAMHandler(args.control_bam, caching=0)
# Regions of Interest
regions = pyDNase.GenomicIntervalSet(args.bedsites)
# Output
treatment_output = open(args.treatment_only_output, "w", buffering=1)
control_output   = open(args.control_only_output, "w", buffering=1)

# Determine Number of CPUs to use
if args.processes:
    CPUs = args.processes
else:
    CPUs = mp.cpu_count()
# NOTE: This roughly scales at about 450mb per 300 regions held in memory
max_regions_cached_in_memory = 50 * CPUs
p = mp.Pool(CPUs)

print "Performing differential footprinting..."

for i in progress.bar(regions):
    # Make sure the interval is actually big enough to footprint to begin with
    if len(i) < 120:
        i.startbp -= 60
        i.endbp += 60
    # Make the optional arguments
    fp_args = {'footprint_sizes': args.footprint_sizes, 'FDR_cutoff': args.FDR_cutoff, 'FDR_iterations': args.FDR_iterations}
    # Perform both comparisons - A against B and B against A
    fp = Diffwell(reads2=reads1, min_score=args.FDR_limit, interval=i, reads=reads2, **fp_args)
    fp2 = Diffwell(reads2=reads2, min_score=args.FDR_limit, interval=i, reads=reads1, **fp_args)
    # Push these tasks to the queue
    p.apply_async(fp, callback=write_treat_to_disk)
    p.apply_async(fp2, callback=write_control_to_disk)
    # Hold here while the queue is bigger than the number of reads we're happy to store in memory
    while p._taskqueue.qsize() > max_regions_cached_in_memory:
        pass
