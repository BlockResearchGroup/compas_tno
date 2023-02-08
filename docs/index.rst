Introduction
=============

.. figure:: /_images/objectives.png
    :figclass: figure
    :class: figure-img img-fluid

.. rst-class:: lead

COMPAS TNO is a Python package to find admissible thrust networks in masonry vaulted structures.

The package implements *Thrust Network Optimisation*, or TNO, within the `COMPAS <https://compas.dev/>`_ framework. TNO is a modular multi-objective optimisation framework to find admissible thrust networks in vaulted masonry structures. Thrust Networks represent the internal forces in masonry structures as a connected force network contained within the structural geometry. Based on the safe theorem of limit analysis, a structure is safe if at least one thrust network is found within its envelope. This set of compressive, internal forces corresponds to a lower-bound equilibrium solution.

With TNO, multiple particular equilibrium states can be obtained, including the structure's minimum and maximum horizontal thrusts, its minimum thickness, its vertical and horizontal collapse loads, and the internal stress following support movements. Beyond providing the internal force distribution for these limit states, an indication of the crack patterns is obtained by looking at the points where the network touches intrados and extrados of the vaults.

Exploring network contained within the structural evelope corresponds to solving a constrained nolinear optimisation problem. TNO sets up, run and output the solution networks to these problems.

To use TNO, check the installation guide on the :ref:`Installation <installation>` and :ref:`Getting Started <gettingstarted>` pages to prepare your machine to run TNO.

A tutorial about the package's elements is provided in :ref:`Tutorial <tutorial>`

A series of simple examples are provided in the :ref:`Examples <examples>` section to start using TNO.

Acknowlegments
--------------

TNO has been developed during `Ricardo Maia Avelino <https://ricardoavelino.github.io>`_'s doctoral research at the `Block Research Group <https://block.arch.ethz.ch>`_, ETH Zurich funded by the Swiss National Science Foundation SNSF (project grant n. 178953).

If you use COMPAS TNO to your research, please refer to our :ref:`publications <publications>`.

References
-----------------

[1] R. Maia Avelino, “Thrust Network Optimisation for the Assessment of Vaulted Masonry Structures”, ETH Zurich, 2023. (to appear)

[2] P. Block, “Thrust Network Analysis: Exploring Three-dimensional Equilibrium”, Massachusetts Institute of Technology, 2009.

[3] J. Heyman, “The stone skeleton”, International Journal of Solids and Structures, vol. 2, no. 2, pp. 249--279, Apr. 1966.


Table of Contents
-----------------

.. toctree::
   :maxdepth: 3
   :titlesonly:

   Introduction <self>
   installation
   gettingstarted
   tutorial
   examples
   api
   license
   publications


.. _References:

