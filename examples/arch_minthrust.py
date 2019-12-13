
from compas_tna.diagrams import FormDiagram
from compas_tno.algorithms import optimise_general

from compas_tno.utilities.constraints import circular_heights
from compas_tno.utilities.constraints import create_cracks
from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_arch

from compas_tno.plotters.plotters import plot_form_xz
from numpy import array

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # All: ['SLSQP', 'PSQP', 'CONMIN', 'COBYLA', 'SOLVOPT', 'KSOPT', 'NSGA2', 'ALGENCAN', 'FILTERSD', 'SDPEN', 'ALPSO', 'ALHSO', 'MIDACO']

    for solver in ['SLSQP', 'PSQP', 'ALGENCAN']:

        # Solvers that reach the optimum: 'SLSQP', 'PSQP', 'ALGENCAN'
        # Solvers that don't achieve maximum (for default parameters): 'KSOPT', 'SDPEN', 'ALPSO', 'ALHSO', 'MIDACO' (Don't use gradient information)
        # Solvers that don't work for the problem: 'CONMIN', 'SOLVOPT'
        # Solvers that return a solution that violate constraints 'COBYLA'
        # Solvers that return error Segmentation fault: 'NSGA2', 'FILTERSD'

        solver_ = 'pyOpt-' + solver

        blocks = 20
        thk = 0.20
        form = create_arch(total_nodes=blocks)
        form = circular_heights(form, thk = thk)
        # form = create_cracks(form, dx =[[0.50, 0.60]], dy = [[-0.1, 0.1]], type = ['top'], view = False)
        # overview_forces(form)
        print_opt = False
        exitflag = 0
        cracks = False

        # Initial parameters

        translation = form.attributes['tmax']
        bounds_width = 5.0
        use_bounds = False
        qmax = 1000
        qmin = -1e-6
        indset = None

        # Optimisation

        exitflag = 1
        count = 0

        # while count < 100 and exitflag is not 0:
        # form = _form(form, keep_q=True)
        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver_,
                                            printout=print_opt,
                                            find_inds=True,
                                            tol=0.01,
                                            translation = translation,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='max',
                                            indset=indset,
                                            bmax = True,
                                            summary=print_opt)

        # Check compression and Save

        q = [attr['q'] for u, v, attr in form.edges(True)]
        qmin  = min(array(q))
        count += 1
        if qmin > -0.1 and exitflag == 0:
            print('Optimisation completed - Trial:',count, 't', thk)
            plot_form_xz(form, radius=0.02, cracks=cracks, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True).show()
