.. _solvers:

********************************************************************************
Solvers
********************************************************************************

The base installation will take care of most of the numerical dependencies of ``TNO`` such as SciPy and NumPy. However a few additional solving libraries must be added to the environment.

From the default installation the ``SLSQP`` solver is already available via SciPy. At least two other packages need to be added: `IPOPT`_ and `SDPT3`_.

IPOPT
=====

IPOPT is a Python-package for Interior-Point-OPTimisation. The installation can be done running the following from the terminal:

.. code-block:: bash

    conda install -c conda-forge cyipopt

More information about the solver and solver options are available here:

* https://pypi.org/project/ipopt/

SDPT3
=====

SDPT3 is a solver to perform conical convex optimisation. It is necessary to perform load-path optimisations. This solver is distributed on the MATLAB's package CVX. Both need to be installed: MATLAB and CVX. The installation is broken in three steps:

1) Download and install MATLAB. Version R2020b is recommended.

* https://www.mathworks.com/products/matlab.html

2) Download and install the CVX to MATLAB:

* http://cvxr.com/cvx/download/

3) Install Python API for MATLAB to control it direclty from Python scripts:

* https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html


Non-mandatory solvers
=====================

Additional solvers can be added to ``TNO`` as listed below:

PyOpt
-----

Library with multiple solvers, some adequate for nonlinear optimisation. It has its own version of ``SLSQP`` and also other solvers including ``GA`` and ``DEVO``. To install run:

.. code-block:: bash

    conda install -c mutirri pyopt

More information abobut the library can be found here:

* http://www.pyopt.org/download.html

PyTorch
-------

The gradients and jacobians are done analytically in ``TNO``. We add an option to compute these through auto-differentiation based on the package ``PyTorch``. The following code should install it:


.. code-block:: bash

    conda install pytorch torchvision -c pytorch
