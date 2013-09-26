__author__ = 'Jason Piper'

#Unfortunately, we have to ensure that the user has numpy is installed,
#as pip is bad at installing numpy and scipy at the same time, and just breaks

try:
    import numpy
except ImportError:
    raise ImportError("Due to a quirk with pip, pyDNase requires numpy to be installed before starting setup")
    
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

setup(
    name='pyDNase',
    version="0.1.0",
    description='DNase-seq analysis library',
    long_description=open('README.rst',"rt").read(),
    author='Jason Piper',
    author_email='j.piper@warwick.ac.uk',
    url='http://jpiper.github.io/pyDNase',
    license='GPLv3',
    ext_modules = [Extension("pyDNase/footprinting/fastbinom", ["pyDNase/footprinting/fastbinom.c"])],
    packages= [
        'pyDNase',
        'pyDNase.footprinting',
    ],

    #Uses a custom version of clint that has a time estimator on the progress bar
    dependency_links = ["http://github.com/jpiper/clint/tarball/develop#egg=clint-0.3.0p"],

    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pysam",
        "clint==0.3.0p"
    ],
    
    package_data = {'pyDNase':["data/*"]},
    
    scripts=[
        "pyDNase/scripts/dnase_average_profile.py",
        "pyDNase/scripts/dnase_to_javatreeview.py",
        "pyDNase/scripts/dnase_wig_tracks.py",
        "pyDNase/scripts/wellington_footprints.py",
        "pyDNase/scripts/examples/example_footprint_scores.py"],
    
    test_suite="test",
)
