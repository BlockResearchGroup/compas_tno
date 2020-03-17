********************************************************************************
Getting Started
********************************************************************************

Create an environment
=====================

We highly recommend to install ``compas_tno`` and related packages in a separate conda environment. In this guide, we will create and use an environmenment based on Python 3.7 with the name "tno", referring to "Thrust Network Analysis", but you can use any other name you like (except for "base", which is the name of the root environment of your conda installation).

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

    conda create -n tno -c conda-forge python=3.7 COMPAS shapely

.. raw:: html

    </div>
    <div class="tab-pane" id="replace_python_osx">

.. code-block:: bash

    conda create -n tno -c conda-forge python=3.7 python.app COMPAS shapely

.. raw:: html

    </div>
    </div>
    </div>
    </div>


Activate the environment
========================

The root environment is active by default. Therefore, you should not forget to activate the "tno" environment whenever you want to work with ``compas_tno`` and its related packages.

.. code-block:: bash

    conda activate tno

Installing parent package compas_tna
===================================

One of the requisites of ``compas_tno`` is the installation alongside its parent package ``compas_tna``.

With pip, :mod:`compas_tna` can be installed directly from the GitHub repo.

.. code-block:: bash

    pip install git+https://github.com/BlockResearchGroup/compas_tna.git#egg=compas_tna


Or from local source files, you need to navigate to your code folder (folder in which you place you packages) and then do:

.. code-block:: bash

    git clone https://github.com/BlockResearchGroup/compas_tna.git
    cd compas_tna
    pip install -e .


Installing compas_tno
=====================

Now the installation of ``compas_tno`` is done in a similar fashion.

With pip, :mod:`compas_tno` can be installed directly from the GitHub repo.

.. code-block:: bash

    pip install git+https://github.com/BlockResearchGroup/compas_tno.git#egg=compas_tno


Or from local source files, you need to navigate to your code folder (folder in which you place you packages) and then do:

.. code-block:: bash

    git clone https://github.com/BlockResearchGroup/compas_tno.git
    cd compas_tno
    pip install -e .


Rhino Configuration
===================

Next, let Rhinoceros know that you installed :mod:`compas_tno`  by typinng in the terminal the following code line:

.. code-block:: bash

    python -m compas_rhino.install -p compas compas_rhino compas_tna compas_tno



Solver's dependency
===================

Note that with the installation of ``compas`` and ``compas_tna`` the packages scipy and numpy are already installed. Thie package scipy contains the nonlinear solver SLSQP that will be the main tool of this package. However, to fully use the functionalities of :mod:`compas_tno` you should add to your installation a series of optimisation packages, that are listed here in order of importance.

1) CVX via MATLAB:

Install MATLAB, version R2019b is recommended:

* https://www.mathworks.com/products/matlab.html

Install CVX as a MATLAB package:

* http://cvxr.com/cvx/download/

Install Python API for MATLAB to control MATLAB direclty from Python scripts:

* https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

You will need CVX installed to run convex optimisation in :mod:`compas_tno`.

2) pyOpt:

Library with extensive list of solvers, some adequate for nonlinear optimisation:

* http://www.pyopt.org/download.html

3) IPOPT:

Library with extensive list of solvers for nonlinear optimisation:

* https://pypi.org/project/ipopt/

