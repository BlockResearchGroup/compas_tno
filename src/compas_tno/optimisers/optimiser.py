


__all__ = ['Optimiser']

class Optimiser(object):

    """The ``Optimiser`` sets the parameters of the optimisation.

    Notes
    -----
    An ``Optimiser`` has the following solvers that can be created.

    *   ``from_SciPy`` : Choses from one of the several options of PyOpt solvers
    *   ``from_PyOpt`` : Choses from one of the several options of PyOpt solvers
    *   ``from_IPOPT`` : Choses from one of the several options of PyOpt solvers
    *   ``from_Matlab`` : Choses from one of the several options of PyOpt solvers
    *   ``from_MMA`` : Choses from one of the several options of PyOpt solvers

    The main parameters to be set in the optimiser vary from the solver and generally are are:

        *  'library'           : 'Scipy',
        *  'solver'            : 'slsqp',
        *  'objective'         : 'min-thrust',
        *  'constraints'       : ['funicular', 'envelope', 'reac_bounds', 'partial_reactions', 'cracks'],
        *  'variables'         : ['ind', 'zb', 'all-qs'],
        *  'use_indset'        : True,
        *  'solver_options'    : {},
        *  'exitflag'          : 1,
        *  'fopt'              : None,
        *  'xopt'              : None,


    To run the optimisation and see results one can use:

        *   ``run``    : Method that will run the specified optimiser in the problem
        *  ``results_summary``:

    """

    __module__ = 'compas_tna.optimiser'


    def __init__(self):
        self.data = {
        'library'           : 'Scipy',
        'solver'            : 'slsqp',
        'objective'         : 'min-thrust',
        'constraints'       : ['funicular', 'envelope', 'reac_bounds', 'partial_reactions', 'cracks'],
        'variables'         : ['ind', 'zb', 'all-qs'],
        'use_indset'        : True,
        'solver_options'    : {},
        'exitflag'          : None,
        'fopt'              : None,
        'xopt'              : None,
        }
        self.fobj = None
        self.fconstr = None
        self.x0 = None
        self.xopt = None



    # This class must separate the functions objective functions and etc...
