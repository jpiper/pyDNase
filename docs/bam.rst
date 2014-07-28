.. _bam:

Getting DNase-seq cut data from BAM files
-----------------------------------------

.. TODO: Make it so loading the example bed and bam file is easy!
.. TODO: Proper class referencing

At the heart of the pyDNase package is the ``BAMHandler`` class, which provides an interfact to the cut data in a BAM file corresponding to a DNase-seq dataset. The interface is extremely simple:

    >>> import pyDNase
    >>> reads = pyDNase.BAMHandler("pyDNase/test/data/example.bam")
    >>> reads["chr6,170863500,170863532,+"]
    {'+': [0,0,0,0,1,0,0,1,1,2,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0],
    '-': [0,10,1,0,1,0,4,9,0,1,0,2,1,0,0,0,0,0,3,0,6,3,0,0,0,1,1,1,3,0,3,6]}

As you can see, querying the ``BAMHandler`` object returns a dictionary containing arrays with cut count on the positive reference strand (+), and cuts on the negative reference strand (-). If you wanted to look at the cuts with reference to something on the opposite strand, you can rotate the data 180 degrees by passing a "-" flag,

    >>> reads["chr6,170863500,170863532,-"]
    {'+': [6,3,0,3,1,1,1,0,0,0,3,6,0,3,0,0,0,0,0,1,2,0,1,0,9,4,0,1,0,1,10,0],
    '-': [0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,1,1,0,0,1,0,0,0,0]}

By default, the BAMHandler caches lookups in 1000bp chunks. You can alter this behaviour at instanstiation. The BAMHandler also gives an interface to the Footprint Occupancy Score (FOS).

.. autoclass:: pyDNase.BAMHandler
    :members: __getitem__, __init__, FOS