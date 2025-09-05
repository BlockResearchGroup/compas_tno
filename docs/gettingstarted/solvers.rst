.. _solvers:

********************************************************************************
Solvers
********************************************************************************

The base installation will take care of most of the numerical dependencies of ``TNO`` such as SciPy and NumPy.
However a few additional solving libraries must be added to the environment.

From the default installation the ``SLSQP`` solver is already available via SciPy.
If needed support to `IPOPT`_ needs to be added.

IPOPT
=====

IPOPT is a Python-package for Interior-Point-OPTimisation.
The installation can be done running the following from the terminal:

.. code-block:: bash

    conda install -c conda-forge cyipopt

More information about the solver and solver options are available here:

* https://pypi.org/project/ipopt/


Loadpath solvers
================

The default loadpath solver is `CLARABEL`_ from `CVXPY`_. Details about the solver are available here:

* https://clarabel.org/stable/

If needed support to `MOSEK`_ needs to be added.

MOSEK
-----

MOSEK is a software package to perform convex optimisation.
It is accessed through `CVXPY`_.
To use MOSEK you need to get yourself a license which can be purchased or obtained for free for scholars.
See their website for more information on how to obtain such license and on their installation process see:

* https://www.mosek.com/products/academic-licenses/

Once you have the licensing worked out we will install MOSEK and CVXPY python packages to use it.

The package CVXPY connects to several convex solvers to solve convex problems.
You can learn about what it does through their website.

* https://www.cvxpy.org/

We can install both CVXPY and MOSEK together with pip:

.. code-block:: bash

    pip install cvxpy mosek

Additional solvers
==================

Additional solvers are coming.