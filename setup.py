import sys
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
from numpy.distutils.misc_util import get_numpy_include_dirs


setup(
    name="TAS-Transmission-Atom-Surface",
    version="0.0.1",
    description="TAS - Transmission spectra with Atom-Surface interactions",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    license="BSD3",
    keywords=["van der Waals", "atom surface interactions", "alkali atoms",
              "spectroscopy", "scientific", "physics", "laser spectroscopy",
              "thin atomic vapours"],
#    url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator",
#    download_url="https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/archive/2.0.5.tar.gz",
    author = 'Tom Peyrot and Nikola Sibalic',
    author_email = 'nikolasibalic@physics.org',

    packages=['tas'],

    zip_safe=False,

)
