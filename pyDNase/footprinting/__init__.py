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

try:
    from . import fastbinom as binom
except ImportError:
    from scipy.stats import binom
    print("""Cannot load Cython binomial CDF implementation.
            loaded the SciPy implementation (Note: This is up to 100x slower!)\n
            note: This happens if you start python from the pyDNase project folder without manually
            compiling the fastbinom extension.""")

from itertools import tee
import numpy.random
import numpy as np
import pyDNase
import warnings
class wellington(object):
    def __init__(self, interval, reads,
                 shoulder_sizes=range(35,36), footprint_sizes = range(11,26,2), FDR=0, bonferroni = 0,):
        self.interval = interval
        #The footprint scores are calculated at instantiation.
        self.scores, self.lengths   = self.calculate(reads, shoulder_sizes,footprint_sizes,FDR, bonferroni)
    def footprints(self, withCutoff=-30, merge = 1):
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
            lbound = max(minimapos-minimaphalffplen,0)
            rbound = min(minimapos+minimaphalffplen,len(templogProb))
            ranges.append((lbound,rbound,templogProb.min()))
            templogProb[lbound:rbound] = 1
            tempMLE[lbound:rbound] = 1
        returnSet = pyDNase.GenomicIntervalSet()
        #Merges overlapping ranges (TODO: documentation)
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
            else:
                merged_ranges = ranges
            #Creates reads GenomicIntervalSet and adds the footprints to them
            for i in merged_ranges:
                rstartbp = self.interval.startbp + i[0]
                #We must add one to the end base of the footprint to account for the BED file format
                rendbp   = self.interval.startbp + i[1] + 1
                region = pyDNase.GenomicInterval(self.interval.chromosome, rstartbp, rendbp, strand="+",score=i[2])
                returnSet += region
        return returnSet

    def window(self, iterable, size):
        """
        Takes reads list (iterable) and returns reads list, each of length size, of rolling windows.
        >>> [i for i in window(range(0,12,2), 3)]
        [(0, 2, 4), (2, 4, 6), (4, 6, 8), (6, 8, 10)]
        """
        iters = tee(iterable, size)
        for i in range(1, size):
            for each in iters[i:]:
                next(each, None)
        return zip(*iters)

    def calculate(self,reads,shoulder_sizes=range(35,36),footprint_sizes = range(11,26,2), FDR=0, bonferroni = 0):
        #TODO: write docstring and doctest

        if self.interval.strand is "-":
            warnings.warn("You're footprinting an interval on the reverse strand! "+
                          "You should be sure you know what you're doing as wellington was not designed for this! "+
                          "Ensure that all the intervals you provide to wellington are on the +ve strand for published behaviour",UserWarning)

        cuts = reads[self.interval]
        forwardArray, backwardArray     = cuts["+"].tolist(), cuts["-"].tolist()

        if FDR:
            numpy.random.shuffle(forwardArray)
            numpy.random.shuffle(backwardArray)

        #Let's compute all the possible binding arrays, this really helps when iterating over multiple footprint sizes
        fw_fpscores_dict = {}
        rv_fpscores_dict = {}

        for fp_size in footprint_sizes:
            halffpround = int((fp_size-1)/2)
            fw_fpscores_dict[fp_size] = [0] * halffpround + [sum(i) for i in self.window(forwardArray, fp_size)]
            rv_fpscores_dict[fp_size] = [0] * halffpround + [sum(i) for i in self.window(backwardArray,fp_size)]

       #Empty list of lists for storing the footprint scores
        log_probs       = [[] for i in range(len(forwardArray))]

        if bonferroni:
            bonferroni_factor = np.log(1.0/sum(reads.samfile.lengths))

        #testing multiple background sizes
        for shoulder_size in shoulder_sizes:
            #This computes the background cut sums for the specified shoulder_size for all basepairs
            f_bindingArray = [0] * (shoulder_size - 1) + [sum(i) for i in self.window(forwardArray,shoulder_size)]
            b_bindingArray = [sum(i) for i in self.window(backwardArray,shoulder_size)] + [0] * (shoulder_size - 1)

            for fp_size in footprint_sizes:
                halffpround = int((fp_size-1)/2)
                #This computes the binding cut sums for the specified fp_size for all basepairs
                fw_fpscores = fw_fpscores_dict[fp_size]
                rv_fpscores = rv_fpscores_dict[fp_size]

                for i in range(shoulder_size+halffpround,len(forwardArray)-shoulder_size-halffpround):
                    xForward  = f_bindingArray[i-halffpround-1]
                    nForward  = xForward + fw_fpscores[i]
                    xBackward = b_bindingArray[i+halffpround+1]
                    nBackward = xBackward + rv_fpscores[i]
                    #This requires that there are DNase Cuts present on both strands
                    if xForward and xBackward:
                        #Null hypothesis for probability to randomly hit background
                        p = float(shoulder_size) / (shoulder_size + fp_size)
                        #This stores the P values and the fp_size used to calculate them in reads tuple. in the log_probs[] list
                        score = binom.logsf(int(xForward - 1), int(nForward), p) + binom.logsf(int(xBackward - 1), int(nBackward), p)
                        log_probs[i].append([score,fp_size])

        #This iterates over all the base pairs in the region and creates arrays for the best score and footprint size
        best_probabilities = []
        best_footprintsizes = []
        for i in range(len(forwardArray)):
            if log_probs[i]:
                best_params =  min(log_probs[i])
                #This catches anything which has floated to negative infinity - but it might not be the best way
                best_score = max(-1000,best_params[0])
                if bonferroni:
                    best_probabilities.append(min(0,best_score - bonferroni_factor))
                else:
                    best_probabilities.append(best_score)
                best_footprintsizes.append(best_params[1])
            else:
                best_probabilities.append(0)
                best_footprintsizes.append(0)
        return (np.array(best_probabilities), np.array(best_footprintsizes))


class wellington1D(wellington):
    def calculate(self,reads,shoulder_sizes=range(35,36),footprint_sizes = range(11,26,2), FDR=0, bonferroni = 0):
        #TODO: write docstring and doctest

        #Here we use some precomputed sums to avoid multiple calculations
        cuts = reads[self.interval]
        forwardArray, backwardArray     = cuts["+"], cuts["-"]
        cutArray     = (forwardArray + backwardArray).tolist()
        if FDR:
            numpy.random.shuffle(cutArray)

        #Empty list of lists for storing the footprint scores
        log_probs       = [[] for i in range(len(cutArray))]

        if bonferroni:
            bonferroni_factor = np.log(1/float(sum(reads.samfile.lengths)))

        for shoulder_size in shoulder_sizes:
            #This computes the background cut  sums for the specified shoulder_size for all basepairs
            f_bindingArray = [0] * (shoulder_size - 1) + [sum(i) for i in self.window(cutArray,shoulder_size)]
            b_bindingArray = [sum(i) for i in self.window(cutArray,shoulder_size)] + [0] * (shoulder_size - 1)
            for fp_size in footprint_sizes:
                halffpround = int((fp_size-1)/2)
                #This computes the binding cut sums for the specified fp_size for all basepairs
                fpscores = [0] * halffpround + [sum(i) for i in self.window(cutArray, fp_size)]
                for i in range(shoulder_size+halffpround,len(forwardArray)-shoulder_size-halffpround):
                    xfor  = f_bindingArray[i-halffpround-1]
                    xback = b_bindingArray[i+halffpround+1]
                    x = xfor + xback
                    n = x + fpscores[i]
                    #This requires that there are actually tags present in these windows
                    if x:
                        p = (shoulder_size*2) / float((shoulder_size*2) + fp_size)
                        #This stores the p values and the fp_size used to calculate them in reads tuple. in the log_probs[] list
                        score = binom.logsf(int(x - 1), int(n), p)
                        log_probs[i].append([score,fp_size])
        #This iterates over all the base pairs in the region and creates arrays for the best score and footprint size
        best_probabilities = []
        best_footprintsizes = []
        for i in range(len(cutArray)):
            if log_probs[i]:
                best_params =  min(log_probs[i])
                #This catches anything which has floated to negative infinity - but it might not be the best way
                best_score = max(-1000,best_params[0])
                if bonferroni:
                    best_probabilities.append(min(0,best_score - bonferroni_factor))
                else:
                    best_probabilities.append(best_score)
                best_footprintsizes.append(best_params[1])
            else:
                best_probabilities.append(0)
                best_footprintsizes.append(0)
        return (np.array(best_probabilities), np.array(best_footprintsizes))