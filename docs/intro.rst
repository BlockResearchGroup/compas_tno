********************************************************************************
Introduction
********************************************************************************

.. figure:: /_images/flow.png
    :figclass: figure
    :class: figure-img img-fluid

``compas_tno`` is a Python package that extends compas_tna for solving the specific problem of assessment of masonry structures according to the Safe Theorem of Limit Analysis. From Heyman [1] and this means that, if a path of compressive and equilibrated forces is found within the structural geometry, the structure is safe under the applied loads and this set of compressive, internal forces corresponds to a lower-bound equilibrium solution.

Block [2] introduced Thrust Network Analysis (TNA), in which the force network is described by its vertical projection in plan, named form diagram, and the equilibrium of the horizontal forces is visualised by a reciprocal graph, named force diagram. This theory is now applied to the openn-source Python package :mod:`compas_tna`.

This package provides the necessary infrastructure to frame TNA as an assessment problem. It sets up and run an optimisation problem in order to find specific networks contained within the bounds of a given masonnry structure.

Check the installation guide on :mod:`Getting Started` and the :mod:`tutorial` for start using :mod:`compas_tno`.

Note that the package is not yet public since it is still in development. If you wish to download or contribute please send me an email (maia@arch.ethz.ch).

References
===================

[1] J. Heyman, “The stone skeleton,” Int. J. Solids Struct., vol. 2, no. 2, pp. 249–279, Apr. 1966.

[2] P. Block, “Thrust Network Analysis,” Massachusetts Institute of Technology, 2009.
