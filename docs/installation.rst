.. _installation:

********************************************************************************
Installation
********************************************************************************

Create a conda environment
==========================

The basic prerequisite to install ``compas_tno`` and the `COMPAS <https://compas.dev>`_ packages is `Anaconda <https://www.anaconda.com/products/individual>`_. Install it on you machine and go to the next step.

We highly recommend to install ``compas_tno`` and related packages in a separate conda environment.

In this guide, we will create and use an environmenment based on Python 3.9 or newer with the name ``tno``.

Open your terminal and type the following to create a new environmenment and install the dependencies `COMPAS <https://compas.dev>`_,  `triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ and the COMPAS Standalone `Viewer <https://github.com/compas-dev/compas_view2>`_:


.. code-block:: bash

    conda create -n tno -c conda-forge python=3.9 COMPAS

The base environment is active by default. Activate the ``tno`` environment with:

.. code-block:: bash

    conda activate tno

Installing compas_tno
=====================

The most direct way to work with ``compas_tno`` is by installing it from PyPI. Navigate in the terminal to the folder and type:

.. code-block:: bash

    pip install compas_tno

This install ``compas_tno`` and the base required COMPAS packages.

Standalone viewer
=================

The Standalone viewer in `COMPAS Masonry <https://github.com/compas-dev/compas_masonry>`_ is used currently to display 3D solutions directly from the terminal. The installation can be done through conda:

.. code-block:: bash

    conda install -c conda-forge compas_masonry

To finalise the installation you need to install a few additonal :ref:`Solvers <solvers>` to your environmenment following the additional guide.

Currently, a work-in-progress UI is being develeoped for Rhino 8 and is coming soon.
