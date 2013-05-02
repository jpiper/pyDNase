================================================
pyDNase - a library for analyzing DNAse-seq data
================================================


Introduction
------------

Many people currently analyzing DNase-seq data are using tools designed for ChIP-seq work, but may be inappropriate for DNAse-seq data where we are less interested in the overlaps of sequenced fragments, but the site at which the cut starts.


Features
-----------

pyDNase allows a user to extract DNase-seq cut data from alignments stored in a BAM file, allowing the fast, random access of DNAse-seq data. The interface makes use of an efficient, compressed memory cache, so that data which is read once is stored in memory (this can be disabled, see documentation). We also provide a reference implementation of the Wellington algorithm, allowing the easy application to any DNase-seq dataset. Note: pyDNase forms part of an analysis toolchain, and thus will not perform peak detection, motif analysis etc.

To identify DNase I hypersensitive sites in DNase-seq data, we recommend using  `HOMER <http://biowhat.ucsd.edu/homer/index.html>`_'s ``findPeaks``
with the following settings
``findPeaks -region -size 500 -minDist 50 -o auto -tbp 0``
followed by merging the overlapping regions with
``bedtools sort -i <input.bed> | bedtools merge -i > <output.bed>``
. We find the results are almost exactly the same as the `HOTSPOT <http://www.uwencode.org/proj/hotspot/>`_ methos employed by ENCODE. See the `HOMER <http://biowhat.ucsd.edu/homer/index.html>`_ documentation for detailed information on how to carry out this procedure.


Examples
------------

We have provided example data files for one genomic region in the ``examples/`` folder.

The most basic functionality we offer is querying a DNase-seq dataset for cut information

    >>> import pyDNase
    >>> reads = pyDNase.BAMHandler("examples/example.bam")
    >>> reads["chr6",170863500:170863532]
    {'+': array([0, 0, 0, 0, 1, 0, 0, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0]),
    '-': array([ 0, 10,  1,  0,  1,  0,  4,  9,  0,  1,  0,  2,  1,  0,  0,  0,  0, 0,  3,  0,  6,  3,  0,  0,  0,  1,  1,  1,  3,  0,  3,  6])}

As you can see, querying the BAMHandler object returns a dictionary containing numpy arrays with cut count on the positive reference strand (+), and cuts on the negative reference strand (-). If you wanted to look at the cuts with reference to something on the opposite strand, you can rotate the data 180 degrees by passing a "-" flag,

    >>> reads["chr6",170863500:170863532,"-"]
    {'+': array([ 6,  3,  0,  3,  1,  1,  1,  0,  0,  0,  3,  6,  0,  3,  0,  0,  0, 0,  0,  1,  2,  0,  1,  0,  9,  4,  0,  1,  0,  1, 10,  0]),
    '-': array([0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 1, 0, 0, 0, 0])}

The BAMHandler loads the DNase-seq cuts in 1000bp blocks, and caches the results in a hash table. You can disable this feature by passing the option ``caching = 0`` when initing the BAMHandler. This is useful if you only anticipate needing to access the reads in a region once in an analysis pipeline, or if you have little RAM.

Often, one may be interested in querying this inxformation for large numbers of regions in the genome (all the DHSs, for eample). We provide a basic way to organise BED files in a GenomicIntervalSet object

    >>> regions = pyDNase.GenomicIntervalSet("examples/example.bed")
    >>> # How many regions are in the BED file?
    ...
    >>> print len(regions)
    1
    >>> # Let's see what the first (and only) interval looks like
    ...
    >>> print regions[0]
    chr6  170863142	170863532	0	0.0	+

Iterating/indexing the GenomicIntervalSet object returns GenomicInterval objects, which are sorted by their score. If you update their scores, then the iteration order will change

We can query the BAMHandler object using a GenomicInterval object

    >>> reads[regions[0]]

If you want to footprint, we provide an easy method to do so. One can import the Wellington object, and get the Wellington footprints for an interval given the reads from a specific experiment.

    >>> from pyDNase.footprinting import wellington
    >>> footprinter = wellington(regions[0],reads)
    >>> footprints = footprinter.footprints(withCutoff=-30)
    >>> print type(footprints)
    <class 'pyDNase.GenomicIntervalSet'>
    >>> for i in footprints:
    ...     print i
    ...
    chr6	170863405	170863454	Unnamed10	-166.028415039	+
    chr6	170863263	170863307	Unnamed8	-143.913098948	+
    chr6	170863339	170863382	Unnamed9	-49.5336608554	+

 These can easily be written to a BED file, for example by

    >>> bedfile = open("footprints.bed","w")
    >>> for i in footprints:
    ...     print >> ofile, i
    >>> bedfile.close()

If you want, you can also extract the raw footprint score (and maximum likelihood estimate of footprint length).

    >>> print footprinter.scores
    [  0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
    0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
    0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00 ...
    >>> print footprinter.lengths
    [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 11 13 11 15 13 11 15 25 25 11
    13 15 17 19 21 23 25 25 11 13 15 17 15 21 19 17 15 13 11 21 25 25 25 25 25 ...


So you can write the raw footprinting scores to a WIG file if you want to using something like

    >>> print "fixedStep\tchrom=" + str(footprinter.interval.chromosome) + "\t start="+ str(footprinter.interval.startbp) +"\tstep=1"
    fixedStep	chrom=chr6	 start=170863142	step=1
    >>> for i in footprinter.scores:
    ...     print i
    0.0
    0.0
    0.0

(obviously you'd need to redirect the output to a file using ``>>`` to actually make a WIG file)

In the ''scripts/'' folder we provide a ``Footprint.py`` script to automatically footprint datasets,  export the results as a WIG track of log probabilities, and export footprints called at a range of sensitivities. All that you need in order to access DNase-seq cut information efficiently and perform footprinting with Wellington is explained above (and scripts to perform common tasks are provided in the ``scripts/`` folder), but the source code is fully documented and I encourage you to have a poke around. If you need any help, email me at j.piper as the first half of my email, and warwick.ac.uk as the second.

Installation
-------------

Before installing, make sure you have the following:

- `Python 2.7.x <http://python.org>`_ (will not work with Python 3, and untested on <2.7)
- A C compiler (`Clang <http://clang.llvm.org/>`_ or `GCC <http://gcc.gnu.org/>`_)
- `samtools <http://samtools.sourceforge.net/>`_

Then, to install pyDNase, simply: ::

    $ pip install pyDNase

It will try and install the following dependencies for you if you haven't already got them.

- NumPy
- Clint
- Pysam

You can also download the contents of this repository and then run: ::

    $ python setup.py install



Contributions
-------------
I highly encourage contributions! This is my first software development project - send any pull requests this way.


Citation
--------
If you use pyDNase or the Wellington algorithm in your work, please cite the following paper:

[Insert paper reference here!]

License
-------
This work is licensed under the GNU GPLv3 license.. ::

    pyDNase - a library for analyzing DNAse-seq data
    Copyright (C) 2013  Jason Piper

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
