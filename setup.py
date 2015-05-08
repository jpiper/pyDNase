__author__ = 'Jason Piper'

import imp
current_version = imp.load_source('lol', 'pyDNase/_version.py').__version__

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

setup(
    name='pyDNase',
    version=current_version,
    description='DNase-seq analysis library',
    long_description=open('README.rst',"rt").read(),
    author='Jason Piper',
    author_email='j.piper@warwick.ac.uk',
    url='http://jpiper.github.io/pyDNase',
    license='GPLv3',
    ext_modules = [Extension("pyDNase.footprinting.WellingtonC", sources = ["pyDNase/footprinting/WellingtonC.c"], extra_compile_args=['-std=c99'])],
    packages= [
        'pyDNase',
        'pyDNase.footprinting',
    ],

    install_requires=[
        # Not enforcing versions for numpy and matplotlib as they can be a bitch to upgrade
        "numpy", # Tested on >=1.5.0
        "matplotlib", # Tested on >=1.2
        "pysam >= 0.7.5",
        "clint >= 0.3.2",
    ],
    
    package_data = {'pyDNase':["data/*"]},
    
    scripts=[
        "pyDNase/scripts/dnase_cut_counter.py",
        "pyDNase/scripts/dnase_average_profile.py",
        "pyDNase/scripts/dnase_to_javatreeview.py",
        "pyDNase/scripts/dnase_wig_tracks.py",
        "pyDNase/scripts/wellington_footprints.py",
        "pyDNase/scripts/wellington_bootstrap.py",
        "pyDNase/scripts/dnase_to_JSON.py",
        "pyDNase/scripts/dnase_ddhs_scorer.py",
        "pyDNase/scripts/examples/example_footprint_scores.py",
        "pyDNase/scripts/dnase_to_JSON.py",
        "pyDNase/scripts/dnase_bias_estimator.py"],
    test_suite="test",
)
