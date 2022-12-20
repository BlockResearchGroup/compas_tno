.. _installation:

********************************************************************************
Installation
********************************************************************************

Create a conda environment
==========================

The basic prerequisite to install ``compas_tno`` and the `COMPAS <https://compas.dev>`_ packages is `Anaconda <https://www.anaconda.com/products/individual>`_. Install it on you machine and go to the next step.

We highly recommend to install ``compas_tno`` and related packages in a separate conda environment.

In this guide, we will create and use an environmenment based on Python 3.8 with the name ``tno``.

Open your terminal and type the following to create a new environmenment and install the dependencies `COMPAS <https://compas.dev>`_,  `triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ and the COMPAS Standalone `Viewer <https://github.com/compas-dev/compas_view2>`_:

.. raw:: html

    <div class="card">
        <div class="card-header">
            <ul class="nav nav-tabs card-header-tabs">
                <li class="nav-item">
                    <a class="nav-link active" data-toggle="tab" href="#replace_python_windows">Windows</a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-toggle="tab" href="#replace_python_osx">OSX</a>
                </li>
            </ul>
        </div>
        <div class="card-body">
            <div class="tab-content">

.. raw:: html

    <div class="tab-pane active" id="replace_python_windows">

.. code-block:: bash

    conda create -n tno -c conda-forge python COMPAS triangle compas_view2

.. raw:: html

    </div>
    <div class="tab-pane" id="replace_python_osx">

.. code-block:: bash

    conda create -n tno -c conda-forge python python.app COMPAS triangle compas_view2

.. raw:: html

    </div>
    </div>
    </div>
    </div>

The base environment is active by default. Activate the ``tno`` environment with:

.. code-block:: bash

    conda activate tno

Installing compas_tno
=====================

The most direct way to work with ``compas_tno`` is by cloning the repository to your local drive. Create a folder to store the code. Navigate in the terminal to the folder and type:

.. code-block:: bash

    git clone https://github.com/BlockResearchGroup/compas_tno.git
    cd compas_tno
    pip install -e .

This install ``compas_tno`` and the base required COMPAS packages.

Standalone viewer
=================

The Standalone viewer `COMPAS View 2 <https://github.com/compas-dev/compas_view2.git>`_ is used currently to display 3D solutions directly from the terminal. The installation can be done through conda:

.. code-block:: bash

    conda install -c conda-forge compas_view2

To finalise the installation you need to install a few additonal :ref:`Solvers <solvers>` to your environmenment following the additional guide.

Currently, a work-in-progress UI is being develeoped for :ref:`Rhino  <rhino>` 6+ and an installation guide is provided.
