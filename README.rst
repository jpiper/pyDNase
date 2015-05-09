================================================
pyDNase - a library for analyzing DNase-seq data
================================================


.. image:: https://travis-ci.org/jpiper/pyDNase.png?branch=master
    :target: https://travis-ci.org/jpiper/pyDNase
.. image:: https://coveralls.io/repos/jpiper/pyDNase/badge.png?branch=master
    :target: https://coveralls.io/r/jpiper/pyDNase?branch=master

Introduction
------------

Many people currently analyzing DNase-seq data are using tools designed for ChIP-seq work, but may be inappropriate for DNase-seq data where one is less interested in the overlaps of sequenced fragments, but the site at which the cut occurs (the 5' most end of the aligned sequence fragment).

We have developed `pyDNase` to interface with a sorted and indexed BAM file from a DNase-seq experiment, allowing efficient and easy random access of DNase-seq cut data from any genomic location, e.g.

    >>> import pyDNase
    >>> reads = pyDNase.BAMHandler(pyDNase.example_reads())
    >>> reads["chr6,170863500,170863532,+"]
    {'+': [0,0,0,1,0,0,1,1,2,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1],
     '-': [0,10,1,0,1,0,4,9,0,1,0,2,1,0,0,0,0,0,3,0,6,3,0,0,0,1,1,1,3,0,3,6]}

Querying the ``BAMHandler`` object returns a dictionary containing lists with DNase cut counts on the positive reference strand (+), and cuts on the negative reference strand (-). `pyDNase` efficiently caches the cut data queried, so that multiple requests from the same genomic locations do not require repeated lookups from the BAM file (this can be disabled).

`pyDNase` comes with several analysis scripts covering several common use cases of DNase-seq analysis, and also an implementation of the Wellington and Wellington 1D footprinting algorithms.

to install `pyDNase`, run::

    $ pip install pyDNase

for full documentation go to: http://pythonhosted.org/pyDNase/


Support
-------

If you're having any troubles, please send an email to `j.piper@warwick.ac.uk` and I'll do my best to help you out. If you notice any bugs, then please raise an issue over at the github repo.


Contributions
-------------
I highly encourage contributions! This is my first software development project - send any pull requests this way. I'm particularly interested in cool analysis scripts that anyone has written.

Reference
---------

.. note ::
    If you use `pyDNase` or the Wellington algorithm in your work, please cite the following paper.
    
    Piper et al. 2013. *Wellington: A novel method for the accurate identification of digital genomic footprints from DNase-seq data*, Nucleic Acids Research 2013; doi: 10.1093/nar/gkt850

License
-------

Copyright (C) 2015 Jason Piper. This work is licensed under the GNU GPLv3 license, see ``LICENCE.TXT`` for details.

