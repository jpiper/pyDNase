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

from . import _version
__version__ = _version.__version__

import os, math, pysam
from clint.textui import progress, puts_err
import sqlite3 as lite
import tempfile, warnings, pickle

def example_reads():
    """
    returns the path to the example BAM file
    """
    return os.path.join(os.path.join(os.path.dirname(__file__), "data"),"example.bam")

def example_regions():
    """
    returns the path to the example BED file
    """
    return os.path.join(os.path.join(os.path.dirname(__file__), "data"),"example.bed")


class BAMHandler(object):
    """
    The object that provides the interface to DNase-seq data help in a BAM file
    """
    def __init__(self, filePath,caching=True,chunkSize=1000):
        """Initializes the BAMHandler with a BAM file

        Args:
            filePath (str): the path of a sorted, indexed BAM file from a DNase-seq experiment
        Kwargs:
            chunkSize (int): and int of the size of the regions to load if caching (default: 1000)
            caching (bool): enables or disables read caching (default: True)
        Raises:
            IOError
        """
        try:
            self.samfile = pysam.Samfile(filePath)
        except IOError:
            errorString = "Unable to load BAM file:{0}".format(filePath)
            raise IOError(errorString)

        #Initialise the empty DNase cut cache with the chromosome names from the BAM file
        self.purge_cache()
        #Do not change the CHUNK_SIZE after object instantiation!
        self.CHUNK_SIZE  = chunkSize
        self.CACHING     = caching

    def __addCutsToCache(self,chrom,start,end):
        """Loads reads from the BAM file into the cutCache. Will not check if reads are already there.

        Args:
            chrom (str): The chromosome
            start (int): The start of the interval
            end (int): The end of the interval
        """
        for alignedread in self.samfile.fetch(chrom, max(start, 0), end):
            if alignedread.is_reverse:
                a = int(alignedread.aend)
                if a <= end +1:
                    self.cutCache[chrom]["-"][a] = self.cutCache[chrom]["-"].get(a, 0) + 1
            else:
                a = int(alignedread.pos) -1
                if a >= start:
                    self.cutCache[chrom]["+"][a] = self.cutCache[chrom]["+"].get(a, 0) + 1
        self.lookupCache[chrom].append(start)

    def __lookupReadsUsingCache(self,startbp,endbp,chrom):
        """Looks up the DNase cut information from the cutCache and returns as a dictionary (private method)

        Args:
            startbp (int): The start of the interval
            endbp (int): The end of the interval
            chrom (str): The chromosome

        """
        #Expand the region to the nearest CHUNK_SIZE and load these reads if they aren't in the cache
        lbound = int(math.floor(startbp / float(self.CHUNK_SIZE)) * float(self.CHUNK_SIZE))
        ubound = int(math.ceil(endbp / float(self.CHUNK_SIZE)) * float(self.CHUNK_SIZE))
        for i in range(lbound,ubound,self.CHUNK_SIZE):
            if i not in self.lookupCache[chrom]:
                self.__addCutsToCache(chrom,i,i+self.CHUNK_SIZE)
        #Fills in with zeroes where the hash table contains no information for each strand.
        fwCutArray  = [self.cutCache[chrom]["+"].get(i, 0) for i in range(startbp,endbp)]
        revCutArray = [self.cutCache[chrom]["-"].get(i, 0) for i in range(startbp,endbp)]
        return {"+":fwCutArray,"-":revCutArray}

    def __lookupReadsWithoutCache(self,startbp,endbp,chrom):
        """Loads reads from the BAM file directly from disk, ignoring the cache (private method)

        Args:
            startbp (int): The start of the interval
            endbp (int): The end of the interval
            chrom (str): The chromosome
        """
        tempcutf = {}
        tempcutr = {}
        for alignedread in self.samfile.fetch(chrom, max(startbp, 0), endbp):
            if alignedread.is_reverse:
                a = int(alignedread.aend)
                if a <= endbp +1:
                    tempcutr[a] = tempcutr.get(a, 0) + 1
            else:
                a = int(alignedread.pos) - 1
                if a >= startbp:
                    tempcutf[a] =tempcutf.get(a, 0) + 1
        fwCutArray  = [tempcutf.get(i, 0) for i in range(startbp,endbp)]
        revCutArray = [tempcutr.get(i, 0) for i in range(startbp,endbp)]
        return {"+":fwCutArray,"-":revCutArray}

    def purge_cache(self):
        """
        Wipes the internal cache - useful if you need finer grain control over caching.
        """

        #Initialise the empty DNase cut cache with the chromosome names from the BAM file
        self.cutCache = {}
        #This helps us differentiate between what's been looked up and when there's just no reads
        self.lookupCache = {}
        for i in self.samfile.references:
            self.cutCache[i]    = {"+":{},"-":{}}
            self.lookupCache[i] = []

    def get_cut_values(self,vals):
        """Return a dictionary with the cut counts. Can be used in two different ways:

        You can either use a string or a GenomicInterval to query for cuts.
        Returns reads dict with "+" corresponding to the +ve strand and "-" has the data with the -ve strand (rotated 180 degrees)

        Args:
            vals: either a string with the format "chr18,500:600,+" or a GenomicInterval object

        >>> BAMHandler(example_reads())["chr6,170863142,170863150,+"]
        {'+': array([ 1,  0,  0,  0,  1, 11,  1,  0]), '-': array([0, 1, 0, 0, 1, 0, 0, 1])}
        >>> BAMHandler(example_reads())["chr6,170863142,170863150,-"]
        {'+': array([1, 0, 0, 1, 0, 0, 1, 0]), '-': array([ 0,  1, 11,  1,  0,  0,  0,  1])}
        """

        if isinstance(vals, GenomicInterval):
            chrom   = vals.chromosome
            startbp = vals.startbp
            endbp   = vals.endbp
            flip    = vals.strand

        elif isinstance(vals, str):
            try:
                chrom,startbp,endbp,flip   = vals.split(",")
                startbp = int(startbp)
                endbp   = int(endbp)
                assert(flip in ["+", "-"])
            except:
                raise ValueError("Malformed query string")

        else:
            raise TypeError("Lookup must be a string or a GenomicInterval")

        if self.CACHING:
            retval = self.__lookupReadsUsingCache(startbp,endbp,chrom)
        else:
            retval = self.__lookupReadsWithoutCache(startbp,endbp,chrom)

        if flip is "-":
            retval["+"], retval["-"] = retval["-"][::-1], retval["+"][::-1]
        return retval

    def __getitem__(self,vals):
        """
        Wrapper for get_cut_values
        """
        return self.get_cut_values(vals)

    def FOS(self,interval,bgsize=35):
        """Calculates the Footprint Occupancy Score (FOS) for a Genomicinterval. See Neph et al. 2012 (Nature) for full details.
        
        Args:
            interval (GenomicInterval): The interval that you want the FOS for

        Kwargs:
            bgsize (int): The size of the flanking region to use when calculating the FOS (default: 35)

        Returns:
            A float with the FOS - returns -1 if it can't calculate it
        """

        cuts = self["{0},{1},{2},{3}".format(interval.chromosome,interval.startbp-bgsize,interval.endbp+bgsize,interval.strand)]
        forwardArray, backwardArray     = cuts["+"], cuts["-"]
        cutArray     = [forwardArray[i] + backwardArray[i] for i in range(len(forwardArray))]

        leftReads   = float(sum(cutArray[:bgsize]))
        centreReads = float(sum(cutArray[bgsize:-bgsize]))
        rightReads  = float(sum(cutArray[-bgsize:]))

        #Here we normalise by region length
        leftReads   /= (bgsize * 1.0)
        centreReads /= (len(interval) * 1.0)
        rightReads  /= (bgsize * 1.0)

        try:
            return ( (centreReads+1.0) / (leftReads + 1.0) ) + ( (centreReads+1.0)/(rightReads + 1.0))
        except BaseException:
            #If it can't calculate the score, return a sentinel value
            return -1


class GenomicIntervalSet(object):
    """Container class which stores and allow manipulations of large numbers of GenomicInterval objects.
    Essentially a way of storing and sorting BED files.
    """

    def __init__(self,filename = None):
        """
        Inits GenomicIntervalSet. You can also specify a BED file path to load the intervals from

        Kwargs:
            filename (str): the path to a BED file to initialize the intervals with

        If no ``filename`` provided, then the set will be empty
        """
        self.intervals = {}

        if filename:
            self.loadBEDFile(filename)

    def loadBEDFile(self,filename):
        """
        Adds all the intervals in a BED file to this GenomicIntervalSet.
        We're quite naughty here and allow some non-standard BED formats (along with the official one):

        chrom chromStart chromEnd
        chrom chromStart chromEnd strand
        chrom chromStart chromEnd name score strand

        Any whitespace (tabs or spaces) will be considered separators, so spaces in names cause a problem!

        .. note::
            If you don't supply a strand, we infer that it's +ve.

        Args:
            filename: the path to a BED file to load

        Raises:
            IOError
        """
        try:
            BEDfile = open(filename, 'rU')
        except IOError:
            errorString = "Cannot load BED file: {0}".format(filename)
            raise IOError(errorString)

        puts_err("Reading BED File...")

        #This is done so that if a malformed BED record is detected, no intervals are loaded.
        records = []
        
        intervalCount = max(enumerate(open(filename)))[0] + 1
        for _ in progress.bar(range(intervalCount)):
            line    = BEDfile.readline()
            #Skip lines in the bed files which are UCSC track metadata or comments
            if not self.__isBEDHeader(line):
                records.append(self.__parseBEDString(line))

        for i in records:
            self.__addInterval(GenomicInterval(i[0], i[1], i[2], i[3], i[4], i[5]))

        BEDfile.close()

    def __malformedBEDline(self,BEDString):
        """
        Raises an exception and prints the offending BED string

        Raises:
            Exception
        """
        #TODO: Make a new exception class, something like malformedBEDException?
        exceptionString = "Malformed BED line: {0}".format(BEDString)
        raise Exception(exceptionString)

    def __isBEDHeader(self,string):
        """
        Returns True/False whether a line in a bed file should be ignored according to
        http://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK
        """
        if string[0] == "#":
            return True
            
        headers = ["name","description","type","visibility","color","itemRgb","useScore","group",
                   "priority","db","offset","maxItems","url","htmlUrl","bigDataUrl","track","browser"]
                   
        for each in headers:
            if string.startswith(each):
                return True
        return False

    def __parseBEDString(self,BEDString):
        """
        Parses the following BED formats
        We're quite naughty here and allow some non-standard BED formats (along with the official one):

        chrom chromStart chromEnd
        chrom chromStart chromEnd strand
        chrom chromStart chromEnd name score strand

        Returns:
            (chrom,startbp,endbp,label,score,strand)

        Raises:
            Exception
        """
        BEDSplit = BEDString.split()

        #Sanity check
        if len(BEDSplit) not in [3,4,6]:
            self.__malformedBEDline(BEDString)

        #Default if only Chrom Start End is detected
        try:
            chrom   = BEDSplit[0]
            startbp = int(BEDSplit[1])
            endbp   = int(BEDSplit[2])
        except:
            self.__malformedBEDline(BEDString)

        label = 0
        score = 0
        strand = "+"

        if len(BEDSplit) is 4:
            if BEDSplit[3] in ["+", "-"]:
                strand = BEDSplit[3]
            else:
                self.__malformedBEDline(BEDString)

        if len(BEDSplit) is 6:
            label  = BEDSplit[3]
            try:
                score = float(BEDSplit[4])
            except ValueError:
                self.__malformedBEDline(BEDString)
            if BEDSplit[5] in ["+", "-"]:
                strand = BEDSplit[5]
            else:
                self.__malformedBEDline(BEDString)

        return chrom,startbp,endbp,label,score,strand

    def __len__(self):
        """
        Return the number of intervals in the set
        """
        intervals = 0
        for each in self.intervals.values():
            intervals += len(each)
        return intervals

    def __iter__(self):
        """
        Iterates over the intervals in the order that the intervals were generated
        """
        for each in sorted(sum(self.intervals.values(),[]), key=lambda peak: peak.importorder):
            yield each

    def __getitem__(self, i):
        """
        Indexes the intervals in the order that the intervals were generated
        """
        return sorted(sum(self.intervals.values(),[]), key=lambda peak: peak.importorder)[i]

    def __delitem__(self,i):
        """
        Deletes an interval from the set using the position i
        """
        pos = self.intervals[self[i].chromosome].index(self[i])
        del self.intervals[self[i].chromosome][pos]

    def __iadd__(self, other):
        """
        Adds all the intervals from an other GenomicIntervalSet or GenomicInterval to this one.

        Args:
            other: either a GenomicInterval or GenomicIntervalSet to be added
        Raises:
            TypeError: A GenomicInterval or GenomicIntervalSet wasn't supplied.
        """
        if isinstance(other,GenomicIntervalSet):
            for i in other:
                self.__addInterval(i)
        elif isinstance(other,GenomicInterval):
            self.__addInterval(other)
        else:
            raise TypeError("Can only add GenomicInterval or GenomicIntervalSet objects to existing GenomicIntervalSet")
        return self

    def __addInterval(self, other):
        """
        Adds a GenomicInterval to the set

        Args:
            other (GenomicInterval): The GenomicInterval to be added.
        """
        if other.chromosome not in self.intervals:
            self.intervals[other.chromosome] = []
        self.intervals[other.chromosome].append(other)

    def resizeRegions(self,toSize):
        """
        Resized all GenomicIntervals to a specific size

        Args:
            toSize: an int of the size to resize all intervals to
        """
        if not type(toSize) == int:
            ValueError("Can only resize intervals to integers")

        for i in self:
            xamount = toSize-(i.endbp-i.startbp)//2
            i.startbp -= xamount
            i.endbp   += xamount
            if (i.endbp-i.startbp) > toSize*2:
                i.endbp -= 1
    def __str__(self):
        return ''.join(str(i) +"\n" for i in self)

class GenomicInterval(object):
    """
    Basic Object which describes reads region of the genome
    """

    #This counts how many GenomicInterval objects have been created
    counter = 0

    def __init__(self, chrom, start, stop, label = 0,score = 0,strand="+"):
        """
        Initialization routine

        Args:
            chrom (str): the chromosome
            
            start (int): the start of the interval
            
            stop (int): the end of the interval

        Kwargs:
            label: The name of the interval (will be given an automatic name if none entered)
            
            score (float): the score of the interval (default: 0)
            
            strand (str): the strand the interval is on (default: "+")
        """

        self.__class__.counter += 1
        self.importorder = self.__class__.counter

        self.chromosome = str(chrom)
        self.startbp    = int(start)
        self.endbp      = int(stop)
        self.strand     = str(strand)
        self.score      = float(score)

        if self.startbp > self.endbp:
            raise Exception("Start location of GenomicInterval is larger than end location!")

        # This is from reads bygone era where we ordered the intervals by import order
        # self.score      = self.__class__.counter

        #This makes up reads fake name if one doesn't exist in the BED file
        if label:
            self.label = str(label)
        else:
            self.label = "Unnamed{0}".format(self.__class__.counter)

        #This contains anything else you want to store about the interval
        self.metadata = {}

    def __str__(self):
        """
        BED representation of the interval
        """
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.chromosome, self.startbp, self.endbp, self.label, self.score, self.strand)

    def __len__(self):
        """
        Returns the length of the GenomicInterval in basepairs
        """
        return self.endbp - self.startbp

class FASTAHandler(object):
    def __init__(self, fasta_file, vcf_file = None):
        self.ffile = pysam.Fastafile(fasta_file)
        self.conn = None
        if vcf_file:
            self.conn = lite.connect(tempfile.NamedTemporaryFile().name)
            with open(vcf_file, 'r') as f:
                records = [(x[0],x[1],x[3],x[4]) for x in (x.split() for x in f if x[0] != "#")]
                with self.conn:
                    cur = self.conn.cursor()
                    cur.execute("CREATE TABLE SNPS(chr TEXT,pos INT, ref TEXT, mut TEXT)")
                    cur.executemany("INSERT INTO SNPS VALUES(?,?,?,?)",records)
                #Manually remove these, as they potentially could be large in memory
                del(records)
    def sequence(self,interval):
        sequence_string = self.ffile.fetch(interval.chromosome,interval.startbp,interval.endbp).upper()
        if not self.conn:
          #  if interval.strand != "+":
           #     sequence_string = sequence_string[::-1]
            return str(sequence_string)
        else:
            query_string = "SELECT chr, pos - ? - 1 as offset,ref,mut FROM SNPS WHERE chr=? and pos BETWEEN ? and ?"
            snps = self.conn.cursor().execute(query_string,(interval.startbp,interval.chromosome,interval.startbp,interval.endbp)).fetchall()
            sequence_list = [i for i in sequence_string]
            for i in snps:
                if sequence_list[i[1]] != i[2]:
                    warnings.warn("MISMATCH IN REF TO SNP - WHAT HAVE YOU DONE?")
                else:
                    #str needed as sqlite returns unicode
                    sequence_list[i[1]] = str(i[3])
            return "".join(sequence_list)

class BiasCalculator(object):
    def __init__(self,bias_file=None):
        if bias_file is None:
            #Load the genomic IMR90 bias from Shirley
            bias_file = open(os.path.join(os.path.join(os.path.dirname(__file__), "data"),"IMR90_6mer.pickle"))
        self.biasdict = pickle.load(bias_file)
    def bias(self,sequence):
        """
        NOTE:   Because bias is calculated from  the centre of a 6-mer,
                the data will be missing the values for the first and last 3 bases
        """
        #Split sequence into consituent 6-mers
        sequence_chunks = [sequence[i:i+6] for i in range(len(sequence)-5)]
        #Look up values of these 6-mers in the precomputed bias database
        fw_bias =  [float(self.biasdict[i]["forward"])for i in sequence_chunks]
        rv_bias =  [float(self.biasdict[i]["reverse"])for i in sequence_chunks]
        #FIXME: Pickled data should use "+" and "-" and not forward and reverse to prevent confusion here
        #FIXME: Fix the pickled data - the reverse positions are off by one!
        return {"+": fw_bias, "-":rv_bias}

class BAMHandlerWithBias(BAMHandler):
    def __init__(self, sequence_object, *args, **kwargs):
        super(BAMHandlerWithBias, self).__init__(*args, **kwargs)
        self.sequence_data = sequence_object
        self.bias_data     = BiasCalculator()

    def __getitem__(self,interval):
        if not isinstance(interval,GenomicInterval):
            raise TypeError("Sorry, but we only support GenomicInterval querying for the Bias Handler at the moment")
        #Note: This is pretty Hacky!
        interval.startbp -= 3
        interval.endbp   += 3
        #Get the sequence data
        bias_values = self.bias_data.bias(self.sequence_data.sequence(interval))
        interval.startbp += 3
        interval.endbp   -= 3
        bias_values["+"] = bias_values["+"][1:]
        bias_values["-"] = bias_values["-"][:-1]
        cuts = self.get_cut_values(interval)

        #Nomenclature used below is that in the He. et al Nature Methods Paper
        #These are N_i^s - note we are using an entire hypersensitive site, and not 50bp like the paper
        N_i = {"+":sum(cuts["+"]) ,"-":sum(cuts["-"])}

        #bias_values are y_i
        for dir in ("-","+"):
            bias_values[dir] = [float(i)/sum(bias_values[dir]) for i in bias_values[dir]]

        #Stupid pass-by-reference
        predicted_cuts = {"+":cuts["+"][:],"-":cuts["-"][:]}
        #Calculate the predicted counts (nhat_i, which is N_i * y_i) based on n-mer cutting preference
        for dir in ("-","+"):
            #For each base
            for num, val in enumerate(predicted_cuts[dir]):
                #Multiply the total number of observed cuts by the bias value
                predicted_cuts[dir][num] = bias_values[dir][num] * N_i[dir]

        #Now we normalised the observed cuts by the expected
        for dir in ("-","+"):
            #For each base
            for num, val in enumerate(predicted_cuts[dir]):
                #Divide the number of observed cuts by the bias value
                pass
                #predicted_cuts[dir][num] = (cuts[dir][num] + 1.0)  / (val + 1.0)
        if interval.strand == "-":
            # That's numberwang, let's rotate the board!
            predicted_cuts["+"], predicted_cuts["-"] = predicted_cuts["-"][::-1], predicted_cuts["+"][::-1]
        return predicted_cuts
