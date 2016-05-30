# Copyright (C) 2016 Jason Piper - j.piper@me.com
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

import warnings, random
import numpy as np
import pyDNase
from . import WellingtonC

class wellington(object):

    def __init__(self,interval,reads,shoulder_sizes=range(35,36), footprint_sizes = range(11,26,2),FDR_cutoff=0.01,FDR_iterations=100,bonferroni=None):

        #Set up the model parameters
        self.shoulder_sizes  = shoulder_sizes
        self.footprint_sizes = footprint_sizes
        self.FDR_iterations  = FDR_iterations
        self.FDR_cutoff      = FDR_cutoff
        if bonferroni:
            self.bonferroni_factor = np.log(1.0/sum(reads.samfile.lengths))
        else:
            self.bonferroni_factor = None

        #Here we check if the interval is big enough for the parameters give
        if len(interval) < (max(shoulder_sizes)*2) + max(footprint_sizes):
            raise ValueError("The interval you're trying to footprint is smaller than the parameters you've passed.")

        self.interval        = interval
        if self.interval.strand is "-":
            warnings.warn("You're footprinting an interval on the reverse strand! "+
                          "You should be sure you know what you're doing as wellington was not designed for this! "+
                          "Ensure that all the intervals you provide to wellington are on the +ve strand for published behaviour",UserWarning)

        #We don't store the read object anymore as it isn't pickleable due it being a file handle.
        self.reads        = reads[self.interval]

        #"None" indicates that these haven't been calculated yet.
        self.storedscore   = None
        self.storedlengths = None
        self.FDR_Result    = None

    #This is a pretty hacky way to allow this object to be called from a multithreaded standpoint.
    #If you call the instanced object itself, it calculates the scores and the null scores, and returns
    #the entire object.
    def __call__(self):
        self.storedscore, self.storedlengths   = self.calculate()
        self.FDR_Result = self.FDRscore()
        return self

    #So I've made FDR_Value, lengths, and scores properties which are calculated ad hoc - they are no longer calculated
    #at instantiation
    @property
    def FDR_value(self):
        if self.FDR_Result != None:
            return self.FDR_Result
        else:
            self.FDR_Result = self.FDRscore()
            return self.FDR_Result
    @property
    def lengths(self):
        if self.storedlengths != None: return self.storedlengths
        else:
            self.storedscore, self.storedlengths   = self.calculate()
            return self.storedlengths

    @property
    def scores(self):
        if self.storedlengths != None: return self.storedscore
        else:
            self.storedscore, self.storedlengths   = self.calculate()
            return self.storedscore

    def FDRscore(self):
        return WellingtonC.percentile(sum([self.calculate(FDR=1)[0] for i in range(self.FDR_iterations)],[]),self.FDR_cutoff)


    def footprints(self, withCutoff=-30, merge = 1):
        """
        This returns reads GenomicIntervalSet with the intervals retrieved below the specific cutoff applied to the selected data
        """
        #This find the positions of all the ranges below the cutoff using reads new method
        ranges = []
        tempMLE, templogProb = np.array(self.lengths), np.array(self.scores)

        #Here we have some different logic for selecting the summits of footprints
        #TODO: Document this part

        while templogProb.min() < withCutoff:
            minimapos = templogProb.argmin()
            minimafplen = tempMLE[minimapos]
            minimaphalffplen = int(minimafplen)/2
            lbound = max(minimapos-(minimaphalffplen),0)
            rbound = min(minimapos+(minimaphalffplen),len(templogProb))
            ranges.append((lbound,rbound,templogProb.min(),minimafplen))
            templogProb[max(lbound-minimafplen,0):min(rbound+minimafplen,len(templogProb))] = 1


        returnSet = pyDNase.GenomicIntervalSet()
        #Merges overlapping ranges (TODO: documentation)
        if ranges:
            # This change here changes the way we merge footprints from the probability trace
            #TODO: Documentation
            if merge:
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
                rstartbp = self.interval.startbp + i[0]
                #We must add one to the end base of the footprint to account for the BED file format
                rendbp   = self.interval.startbp + i[1] + 1
                region = pyDNase.GenomicInterval(self.interval.chromosome, rstartbp, rendbp, strand="+",score=i[2])
                returnSet += region
        return returnSet

    def calculate(self,FDR=0):
        #The passes the read data to the Cython wrapper for Wellington. This yields a roughly 4x speed improvement
        return WellingtonC.calculate(FDR,self.reads["+"],self.reads["-"],self.footprint_sizes,self.shoulder_sizes, self.bonferroni_factor)


class wellington1D(wellington):
    def calculate(self,FDR=0):
        return WellingtonC.calculate1D(FDR,self.reads["+"],self.reads["-"],self.footprint_sizes,self.shoulder_sizes, self.bonferroni_factor)
