********************************************************************************
Getting Started
********************************************************************************

Create an environment
=====================

We highly recommend to install ``compas_tno`` and related packages in a separate conda environment. In this guide, we will create and use an environmenment based on Python 3.7 with the name "tno", referring to "Thrust Network Optimisation", but you can use any other name you like (except for "base", which is the name of the root environment of your conda installation).

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

    conda create -n tno -c conda-forge python=3.7 COMPAS

.. raw:: html

    </div>
    <div class="tab-pane" id="replace_python_osx">

.. code-block:: bash

    conda create -n tno -c conda-forge python=3.7 python.app COMPAS

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
====================================

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

Note that with the installation of ``compas`` and ``compas_tna`` the packages scipy and numpy are already installed.

The package scipy contains the nonlinear solver SLSQP that will be the main tool of this package. Therefore, if you got this far you can already test most of the examples using SLSQP as the solver.

However, to fully use the functionalities of :mod:`compas_tno` you should add to your installation a series of optimisation packages and helpers, that are listed here in order of importance.


PyTorch (Helper)
================

The derivatives of many of the constraints on the optimisation are computed by autodifferentiation based on the package PyTorch. So to be able to use MMA and IPOPT you will need PyTorch installed. To install it please run from your terminal:


.. code-block:: bash

    conda install -c pytorch pytorch

MMA (Solver)
============

The first additional solver MMA already ships with ``compas_tno``, so if all went good on the installation of PyTorch, congratulations, you can already use MMA as an option for the solving.

Learn more about MMA (Method of Moving Asymptotes) by Krister Svanberg:

* https://people.kth.se/~krille/mmagcmma.pdf

IPOPT (Solver)
==============

IPOPT is a well-developed package for nonlinear optimisation. To install just run from your terminal:

.. code-block:: bash

    conda install -c conda-forge cyipopt

More information about the solver and solver options are available here:

* https://pypi.org/project/ipopt/


SDPT3 (Solver)
==============

To perform load-path optimisation the suggestion is to use SDPT3. This solver, however, is only available on the MATLAB package CVX. That means that to perform the optimisaiton you will need to have both: MATLAB and CVX on your computer. The installation of this package is composed of three steps:

1) Download and install MATLAB. Version at least R2019b is recommended.

* https://www.mathworks.com/products/matlab.html

2) Download and install the MATLAB package CVX to perform convex optimisation:

* http://cvxr.com/cvx/download/

3) Install Python API for MATLAB to control MATLAB direclty from Python scripts:

* https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

If you don't need to perform load path optimisation you can skipt this installation.

PyOpt (Library of solvers)
==========================

Library with extensive list of solvers, some adequate for nonlinear optimisation. It has its own version of SLSQP and also other solvers including GA and DEVO. To install run:

.. code-block:: bash

    conda install -c mutirri pyopt

More information abobut the library can be found here:

* http://www.pyopt.org/download.html

