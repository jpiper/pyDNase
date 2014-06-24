.. _faqs:

Frequently Asked Questions
--------------------------

Here are common questions we get. If there are any questions about pyDNase or general DNase-seq analysis, either raise an issue on GitHub or email me on ``j.piper@warwick.ac.uk``


How can I identify hypersensitive sites in DNase-seq data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To identify DNase I hypersensitive sites in DNase-seq data, we recommend using  `HOMER <http://biowhat.ucsd.edu/homer/index.html>`_'s ``findPeaks`` with the  parameters: ``findPeaks -region -size 500 -minDist 50 -o auto -tbp 0``, converting the HOMER peaks to a BED file using ``pos2bed.pl`` and then merging the overlapping regions with::

    $ bedtools sort -i <input.bed> | bedtools merge -i - > <output.bed>

We find the results are almost exactly the same as the `HOTSPOT <http://www.uwencode.org/proj/hotspot/>`_ method employed by ENCODE. See the `HOMER <http://biowhat.ucsd.edu/homer/index.html>`_ documentation for detailed information on how to carry out this procedure.

pyDNase won't install/import or gives weird errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The most common issue here is that you have old versions of the dependencies - namely ``scipy``, ``numpy``, or ``pysam`` installed - try updating these to their latest version.