
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check
from compas_thrust.algorithms.scale import scale_form

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.utilities.constraints import set_height_constraint
from compas_thrust.diagrams.form import overview_forces
from compas_thrust.utilities.symmetry import create_sym2
from compas_thrust.utilities.symmetry import replicate2
from compas_thrust.utilities.symmetry import create_sym
from compas_thrust.utilities.symmetry import replicate
from compas_thrust.utilities.loads import not_sym_load
from compas_thrust.utilities.loads import fill_load

from compas_thrust.diagrams.form import _form
from compas_thrust.diagrams.form import remove_feet

from compas_thrust.diagrams.form import adapt_tna
from compas_thrust.diagrams.form import remove_feet

from compas_thrust.plotters.plotters import plot_form
# from compas_viewers.meshviewer import MeshViewer

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    lps = []

    for i in range(4,5):

        print('\n\n\n\n------------------ Form ',str(i),'\n\n')

        j = 2

        file = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+str(j)+'_0'+str(i)+'_sym.json'
        file_save = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_complete_maxz.json'
        file_complete = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_complete.json'

        form = FormDiagram.from_json(file)
        overview_forces(form)
        form = set_height_constraint(form)
        # z = [form.vertex_coordinates(key)[2] for key in form.vertices()]
        # zmax = max(z)
        # r = zmax/1.0
        # form = scale_form(form, r)
        # form = not_sym_load(form)
        # form = fill_load(form)
        plot_form(form).show()
        overview_forces(form)
        print('Horizontal checks: {0}'.format(horizontal_check(form)))
        check_constraints(form, show=True)
        
        # Initial parameters

        tmax = 5.0 # form.attributes['tmax']
        bounds_width = 5.0
        use_bounds = False
        qmax = 20.0
        indset = None
        nsol = 5
        sol = 0
        sols = []
        forms = []
        tol = 0.001

        # plot_form(form,radius=0.04).show()

        # Optimisation Routine - Many trials and Many solutions

        for k in range(30):
            # form = fill_load(form)
            form = _form(form, keep_q=True)
            print('Optimisation trial {0}.'.format(k))
            if sol == 0 and k % 5 == 0 and k > 0:
                form = _form(form)
                print('Shuffle the form')
            if k > 0 and k % 5 is not 0:
                qi_max  = max(array(q))
                if qi_max > 0.95 * qmax:
                    qmax = qmax * 1.25
                    print('upgrate on qmax')
                else:
                    bounds_width = bounds_width/k
                    use_bounds = False
            fopt, qopt = optimise_single(form, qmax=qmax, solver=None,
                                            polish='slsqp',
                                            population=600,
                                            generations=200,
                                            printout=None,
                                            tol=0.01,
                                            t = tmax,
                                            opt_max=False,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='loadpath',
                                            indset=indset,
                                            buttress=False)
            q = [attr['q'] for u, v, attr in form.edges(True)]
            print('Result iteration {0}: {1}'.format(k,fopt))
            qmin  = min(array(q))
            if qmin > -0.1 and fopt is not None and check_constraints(form, show= False) < 500.0:
                forms.append(deepcopy(form))
                sols.append(fopt)
                print('Optimisation found a result after {0} trials.'.format(k+1))
                sol = sol + 1
                if (sol >= nsol or abs(fopt-sols[len(sols)-2]) < 0.001) and len(sols) > 1:
                    break

        # Printing Results

        print('Solutions for example: ')
        print(sols)
        n = argmin(sols)
        print('Optimum Solution: {0}'.format(n))
        form = forms[n]
        lps.append(form.attributes['loadpath'])
        form.to_json(file_save)
        reactions(form)
        print('Horizontal checks: {0}'.format(horizontal_check(form)))
        overview_forces(form)

        # Replicate and Save Complete

        # form = FormDiagram.from_json(file_complete)
        form_ = replicate(form, file_complete)
        reactions(form_)
        overview_forces(form_)
        check_constraints(form_)
        form_.to_json(file_save)

    print(lps)

