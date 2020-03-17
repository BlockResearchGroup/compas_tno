********************************************************************************
Getting Started
********************************************************************************


Installation
============

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



Main dependency
===================

The main dependency of the program is the parent package :mod:`compas_tna`. To install it follow the instructions of installation on the website:

* https://blockresearchgroup.github.io/compas_tna/tutorial.html
* https://blockresearchgroup.github.io/compas_tna/examples.html
* https://blockresearchgroup.github.io/compas_tna/api.html


Solver's dependency
===================

Note that with the insntallation of :mod:`compas` the packages scipy and numpy are already installed. Thie package scipy contains the nonlinear solver SLSQP
that will be the main tool of this package. However, to fully use the functionalities of :mod:`compas_tno` you should add to your installation a series of optimisation
packages, that are listed here in order of importance.

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

