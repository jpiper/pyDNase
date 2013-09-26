.. _intervals:

Handling Genomic Intervals
---------------------------

``GenomicInterval``
~~~~~~~~~~~~~~~~~~~

The ``GenomicInterval`` is effectively pyDNase's way of storing a BED interval. There are three mandantory fields when creating a new ``GenomicInterval``::

    >>> import pyDNase
    >>> interval = pyDNase.GenomicInterval("chr1",100,200)
    >>> print interval
    chr1	100	200	Unnamed1	0.0	+


.. autoclass:: pyDNase.GenomicInterval
    :members: __init__


You might be wondering why this by itself is helpful. It isn't, until you consider that you can use collections of multiple ``GenomicInterval`` instances in a ``GenomicIntervalSet``

``GenomicIntervalSet``
~~~~~~~~~~~~~~~~~~~~~~

Often, one may be interested in querying cut information for large numbers of regions in the genome (all the DHSs, for example). We provide a basic way to organise BED files using a ``GenomicIntervalSet`` object. 

    >>> import pyDNase
    >>> regions = pyDNase.GenomicIntervalSet("pyDNase/test/data/example.bed")
    >>> print len(regions)  # How many regions are in the BED file?
    1
    >>> print regions
    chr6	170863142	170863532	0	0.0	+

Iterating/indexing the GenomicIntervalSet object returns GenomicInterval objects, which are sorted by their order of creation (so the order of the BED file if importing a BED file). You can sort by any of the other attributes that the GenomicInterval has, for example, to iterate by score,

    >>> for i in sorted(regions,key=lambda x: x.score):
            print i

The key here, is that as well as querying the BAMHandler for cuts using a string, we can also query using a GenomicInterval object

    >>> reads = pyDNase.BAMHandler("pyDNase/test/data/example.bam")
    >>> reads[regions[0]]                                                   #Note: I've truncated this output
    {'+': array([1,0,0,0,1,11,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,2, ...]),
     '-': array([0,1,0,0,1,0 ,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,5,0,0, ...])}
            
For example, one could use this to efficiently calculate the total number of cuts in a DNase-seq dataset using the intervals in a BED file

    >>> readcount = 0
    >>> for interval in regions:
            readcount += reads[interval]["+"].sum() + reads[interval]["-"].sum()
    >>> print readcount
    3119

We have overloaded the ``+`` operator you can directly add other ``GenomicIntervalSet`` or ``GenomicInterval`` objects, and you can delete intervals using the ``del`` keyword thus:

    >>> print regions
    chr6	170863142	170863532	0	0.0	+

    >>> regions += pyDNase.GenomicInterval("chr10","100000000","200000000", "0", 10, "-")
    >>> print regions
    chr6	170863142	170863532	0	0.0	+
    chr10	100000000	200000000	0	10.0    -

    >>> del regions[0]
    >>> print regions
    chr10	100000000	200000000	0	10.0	-

.. autoclass:: pyDNase.GenomicIntervalSet
    :members: __init__, resizeRegions, loadBEDFile


.. _python: http://www.python.org/
.. _samtools: http://samtools.sourceforge.net/
.. _homebrew: http://brew.sh/
.. _cython: http://www.cython.org
.. _NumPy: http://www.numpy.org/‎
.. _clint: https://github.com/kennethreitz/clint
.. _pysam: https://code.google.com/p/pysam/
.. _SciPy: http://www.scipy.org/‎
.. _matplotlib: http://www.matplotlib.org
.. _pip: https://pypi.python.org/pypi/pip