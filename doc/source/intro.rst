Solution method and assumptions
===============================

Problem statement
-----------------

This code computes the coherent transmission through a slab of
atomic vapour contained in thin vapour cell.


Solution method
---------------

Incident field and the field radiated by the atoms interfere
to produce the output field.
Model accounts for the velocity of the atoms by resolving Bloch equation
introducing the hydrodynamic derivative and assuming Maxwell-Boltzmann velocity
distribution. Boundary condition for the limit is a loss of coherence
at surfaces. The atom-surface interaction is regularly accounted for
introducing a long-range atom-surface potential and solving the coherence
equation using the constant variation method.


Assumptions
-----------

The model is valid for low enough absorption since it assume that the driving
field is not depleted along the cell. (First order density derivation of
the transmission - Born Approximation.)
