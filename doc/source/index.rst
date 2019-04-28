TAS - Transmission spectra with Atom-Surface interactions documentation
=======================================================================

A Python package :obj:`tas` calculates transmission spectra through
optically thin slabs of atomic vapours contained in
spectroscopic cells including important effects of  thin cells:

- two closely spaced walls of the spectroscopic cell act
  as low-finesse cavity for driving light and emitted  light;

- atoms colliding with the cell walls depolarise, and their subsequent
  transient dynamics also changes transmission spectra; finally

- **atom-surface van der Waals interaction** induces level shifts.
  For thickness smaller than λ/(2π), where
  λ is the wavelength of the strongest dipole allowed transition
  from each of the two energy levels whose transition spectrum we
  explore, the shift of the transition energy
  is :math:`\propto C_3/R^3`, where :math:`C_3`
  is van der Waals constant for given atomic transition.

Generated  spectra can be used for fitting results
of atom-vapour laser spectroscopy to obtain strength of :math:`C_3` for
atom-surface interactions. Detailed  documentation follows bellow.
For detail on the physics and measurements **please cite
original paper**.


.. toctree::
   :maxdepth: 2

   intro
   installation
   tas
   examples

Credits
--------

:Authors:
    Tom Peyrot and Nikola Šibalić

:Licence: BSD-3

:Version: 0.0.1 of 2019/04/17
