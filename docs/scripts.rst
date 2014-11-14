.. _scripts:

Useful Scripts
--------------

pyDNase installs several scripts which also serve as examples of how to use the pyDNAse API. After you have installed pyDNase these will all be installed into your $PATH, so you can run these directly from your terminal by typing in the name of the script (including the .py extension)

Please have a rummage through the source - it's all documented (and hopefully understandable!)


example_footprint_scores.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script tests that everything has been installed and will run correctly. Upon running it, you should see the following window

.. image:: images/example_footprint_scores_output.png
    
If so, congratulations! Everything has installed properly.
The red and blue bars correspond to cuts on the positive and negative strand,
respectively, and the black line represents the raw Wellington footprint scores.

dnase_cut_counter.py
~~~~~~~~~~~~~~~~~~~~~~~~

This is a very simple script that take a BED and a BAM file and produce a new BED file where the score is the number of DNase cuts found in this region.

.. program-output:: python ../pyDNase/scripts/dnase_cut_counter.py -h

dnase_average_profile.py
~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: images/K562CTCFCHIP.png


Average profile of DNase I activity surrounding ChIP-seq confirmed CTCF sites in K562 data.

Average profile plots illustrating DNase activity surrounding a set of regions are frequently used in papers.
Here, we provide a simple way to generate one

.. program-output:: python ../pyDNase/scripts/dnase_average_profile.py -h

Hopefully this is self-explanatory. This script uses matplotlib to generate the output,
so it will write a filetype based on the file extension provided (e.g. ``out.png`` or ``output.pdf``).

dnase_wig_tracks.py
~~~~~~~~~~~~~~~~~~~~~~

Often, we want to visualise the raw cut data (just the 5' most ends of the cuts) from a DNase-seq experiment, as visualising the pileups isn't helpful here. Here's the FMR1 promoter viewed as a BAM file in IGV

.. image:: images/FMR1a.png

and here's the corresponding cut locations.

.. image:: images/FMR1b.png

We provide ``dnase_wig_tracks.py`` that generates a WIG file (we recommend you convert it to a BigWIG file using UCSC's `wigToBigWig`)
based on a BAM file a list of regions of interest

.. program-output:: python ../pyDNase/scripts/dnase_wig_tracks.py -h

Note that by default, cuts on the reverse strand will be reported as negative numbers (for visualisation). If you want to be using this data for something else, you can pass the ``-r`` flag, which will use the real number of cuts.

dnase_to_JSON.py
~~~~~~~~~~~~~~~~

This outputs DNase-seq data to JSON

.. program-output:: python ../pyDNase/scripts/dnase_to_JSON.py -h


dnase_to_javatreeview.py
~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: images/K562AP1CHIP.png

Want to make a heatmap? Love JavaTreeView_? So do we! This script will generate a CSV file that you can put straight into JavaTreeView to visualize your data.

The options to be aware of here are ``-i`` and ``-a``

.. program-output:: python ../pyDNase/scripts/dnase_to_javatreeview.py -h


wellington_footprints.py
~~~~~~~~~~~~~~~~~~~~~~~~

So you want to get footprints from your data? No problem. We provide a handy script that will do this for you. There's lots of options here, so please read through them carefully. The most basic usage of the script uses the default parameters described in our original paper. If anything goes wrong at any point, then there should be useful error messages telling you exactly what went wrong.

.. program-output:: python ../pyDNase/scripts/wellington_footprints.py -h


.. _JavaTreeView: http://jtreeview.sourceforge.net/

