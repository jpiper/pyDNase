.. _footprinting:

Footprinting DNase-seq data
---------------------------

.. note :: We provide wellington_footprints.py as a script, which will automate footprinting for end users. This below information only necessary if you want to do something fancy. You might want to read the documentation in the source for more information.

We provide a simple interface for footprinting in the ``pyDNase.footprinting`` module. There are two footprinters, ``pyDNase.footprinting.wellington`` and ``pyDNase.footprinting.wellington1D``, which inherits from wellington and overrides the ``calculate`` method with a 1D version.

If you want to footprint, we provide an easy method to do so. One can import the Wellington object, and get the Wellington footprints for an interval given the reads from a specific experiment.

    >>> import pyDNase
    >>> import pyDNAse.footprinting as fp
    >>> regions = pyDNase.GenomicIntervalSet("pyDNase/test/data/example.bed")
    >>> reads = pyDNase.BAMHandler("pyDNase/test/data/example.bam")
    >>> footprinter = fp.wellington(regions[0],reads)
    >>> footprints = footprinter.footprints(withCutoff=-30)
    print footprints
    chr6	170863264	170863306	Unnamed4	-150.07397301	+
    chr6	170863338	170863383	Unnamed5	-47.9227745068	+
    chr6	170863404	170863454	Unnamed6	-164.119817804	+

These can easily be written to a BED file, for example by

    >>> with open("output.bed","w") as bedout:
    >>>     bedout.write(str(footprints))

If you want, you can also extract the raw footprint score.

    >>> print footprinter.scores
    [  0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
    0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
    0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00 ...

So you can write the raw footprinting scores to a WIG file if you want to using something like

    >>> print "fixedStep\tchrom=" + str(footprinter.interval.chromosome) + "\t start="+ str(footprinter.interval.startbp) +"\tstep=1"
    fixedStep    chrom=chr6     start=170863142    step=1
    >>> for i in footprinter.scores:
    ...     print i
    0.0
    0.0
    0.0

(you'd need to redirect these ``print`` statements to a file object to write the actualy WIG file)