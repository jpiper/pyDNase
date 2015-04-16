.. _installation:

Installation
------------

Supported systems
~~~~~~~~~~~~~~~~~

:Hardware:

    1GB of RAM and a 64-bit operating system. RAM usage heavily depends on what you're doing, but 1GB is the bare minimum if you disable caching.

:Software:

   Tested on OS X 10.8 and on Ubuntu 11.10, but should run fine on any other *NIX flavour as long as the prerequisites are fulfilled.


   .. warning::
        Windows is not supported.


Pre-installation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to install :mod:`pyDNase`, the following software is required. Most people will already have most of these on their system. I have attempted to list them in the order that you need to install them in.

#. A compiler suite
    You can check by opening up the terminal and typing ::

        $ clang --version

    or ::

        $ gcc --version

    As long as you get a response from one of these, you're good to go. Failing that...
    
    * **On OS X < 10.7.3:** Install "Xcode" from https://developer.apple.com/downloads/
    * **On OS X >= 10.7.3:** Install "Command Line Tools for Xcode" from https://developer.apple.com/downloads/
    (you can also install Xcode, but this is overkill)
    * **On Ubuntu:** Install with ``sudo apt-get install build-essentials``
   
   .. note::
        If you're using another *NIX distro, I assume you know what you're doing.

#. Python_ >= 2.6 (including Python 3!)
    * This will come installed with OS X or any respectable \*NIX distro.

#. pip_:
        Used for automated installation of Python packages. If you don't already have pip_ installed, you can use the following command to install it ::

            $ curl https://raw.github.com/pypa/pip/master/contrib/get-pip.py | python

#. samtools_
    * **On OS X** the simplest way to install samtools_ is using the homebrew_ command ``brew tap homebrew/science`` followed by ``brew install homebrew/science/samtools``.
    * **On Ubuntu** you can use ``sudo apt-get install samtools``

#. bedtools_
    * **On OS X** the simplest way to install bedtools_ is using the homebrew_ command ``brew tap homebrew/science`` followed by ``brew install homebrew/science/bedtools``.
    * **On Ubuntu** you can use ``sudo apt-get install bedtools``


Optional installs
~~~~~~~~~~~~~~~~~

#. pybedtools_
    * This is only required if you want to use the ``dnase_bias_estimator.py``. Provided you installed pip_, you should be able to simply run ::

        $ pip install pybedtools


Installing :mod:`pyDNase`
~~~~~~~~~~~~~~~~~~~~~~~~~

To install, simply ::

    $ pip install pyDNase

This will attempt to download, compile, and install the python dependencies (``clint``, ``numpy``, ``pysam``, and ``matplotlib``) automatically. However, due to a myriad of reasons it might not work. If this is the case, go and install these manually in said order, then try ``pip install pyDNase`` once more.

.. _pybedtools: https://pythonhosted.org/pybedtools/
.. _python: http://www.python.org/
.. _samtools: http://www.htslib.org/
.. _bedtools: http://bedtools.readthedocs.org/en/latest/
.. _homebrew: http://brew.sh/
.. _NumPy: http://www.numpy.org/‎
.. _clint: https://github.com/kennethreitz/clint
.. _pysam: https://code.google.com/p/pysam/
.. _SciPy: http://www.scipy.org/‎
.. _matplotlib: http://www.matplotlib.org
.. _pip: https://pypi.python.org/pypi/pip
