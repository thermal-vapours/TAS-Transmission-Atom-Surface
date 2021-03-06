TAS - Transmission spectra with Atom-Surface interactions
=========================================================

A Python package for generating transmission spectra through
optically thin slabs of atomic vapours contained in
spectroscopic cells of thickness < λ/(2π), where
atom-surface van der Waals interaction induces level shifts.
Generated  spectra can be used for fitting results
of atom-vapour laser spectroscopy to obtain strength of
atom-surface interactions.

For more details on problem statement, solution method,
assumptions and general use of package
[see documentation](https://tas-transmission-atom-surface.readthedocs.io/en/latest/).

For example of use of generated spectra to obtain
van der Waals interaction coefficient between alkali
atoms and sapphire surface, [see paper](https://doi.org/10.1103/PhysRevA.100.022503).
Example of using data from paper and fitting procedure
to obtain van der Waals interaction coefficient
is provided in [example notebook](/examples/basic_example.ipynb)
that can be [run in browser using Binder](https://mybinder.org/v2/gh/thermal-vapours/TAS-Transmission-Atom-Surface.git/master?urlpath=lab%2Ftree%2Fexamples%2Fbasic_example.ipynb).

[![Documentation Status](https://readthedocs.org/projects/tas-transmission-atom-surface/badge/?version=latest)](https://tas-transmission-atom-surface.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/thermal-vapours/TAS-Transmission-Atom-Surface.svg?branch=master)](https://travis-ci.org/thermal-vapours/TAS-Transmission-Atom-Surface) [![PyPI version](https://badge.fury.io/py/TAS-Transmission-Atom-Surface.svg)](https://badge.fury.io/py/TAS-Transmission-Atom-Surface) [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/thermal-vapours/TAS-Transmission-Atom-Surface.git/master?urlpath=lab%2Ftree%2Fexamples%2Fbasic_example.ipynb)

Authors
-------

Tom Peyrot and Nikola Šibalić, 2019

**Please cite as**: [T. Peyrot, N. Šibalić, Y.R.P. Sortais, A. Browaeys, A. Sargsyan, D. Sarkisyan, I.G. Hughes, C.S. Adams, "Measurement of the atom-surface van der Waals interaction by transmission spectroscopy in a wedged nano-cell", *Phys. Rev. A* **100**, 022503 (2019)](https://doi.org/10.1103/PhysRevA.100.022503)

License
-------

All the files distributed with this program are provided subject to the
BSD-3-Clause license. A copy of the license is provided.
