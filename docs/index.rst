Introduction
=============

.. figure:: /_images/flow.png
    :figclass: figure
    :class: figure-img img-fluid

``compas_tno`` is a Python package used for assessing masonry structures according to the Safe Theorem of Limit Analysis. It works by finding admissible stress states, which are defined as a path of compressive and equilibrated forces contained within the structural geometry. Therefore, from the Safe theorem, if one admissible stress state is found, the structure is safe under the applied loads. This set of compressive, internal forces corresponds to a lower-bound equilibrium solution (:ref:`Heyman, 1966<References>`).

Block (:ref:`2009<References>`) introduced Thrust Network Analysis (TNA) which describes the structure as a network with internal axial forces in the edges and external loads applied to its vertices. Using the principles of graphic statics (:ref:`Wolfe, 1921<References>`), this network is defined by its vertical projection in the plan (the *form diagram*) and the horizontal equilibrium is visualised by a reciprocal graph (the *force diagram*). TNA has been used as an interactive design tool and is implemented in the open-source Python package :mod:`compas_tna`.

This package provides the necessary infrastructure to frame TNA as an assessment problem. It sets up and runs optimisation problems to find specific networks contained within the bounds of a given masonry structure. With these specific networks, one can assess the level of stability of the structure, which goes beyond merely proving that the structure is safe. Additional information about the package content can be found under :mod:`Publications`.

This package is in development by the Block Research Group in a project funded by the Swiss National Science Foundation SNSF (project grant n. 178953). Thus, it has not yet been released to the public. If you wish to download or contribute please an email to (maia@arch.ethz.ch).

Check the installation guide on :mod:`Getting Started` and the :mod:`Tutorial` for start using ``compas_tno``.


Table of Contents
=================

.. toctree::
   :maxdepth: 3
   :titlesonly:

   Introduction <self>
   gettingstarted
   tutorial
   examples
   api
   license
   publications


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

.. _References:

References
===========

[1] P. Block, “Thrust Network Analysis,” Massachusetts Institute of Technology, 2009.

[2] J. Heyman, “The stone skeleton,” Int. J. Solids Struct., vol. 2, no. 2, pp. 249–279, Apr. 1966.

[3] W.S. Wolfe, Graphical analysis: a handbook on graphic statics. New Cork:McGraw-Hill, 1921.
