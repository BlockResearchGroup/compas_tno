
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check

from compas_thrust.utilities.utilities import check_constraints
from compas_thrust.utilities.utilities import oveview_forces
from compas_thrust.utilities.utilities import create_sym_2
from compas_thrust.utilities.utilities import replicate2

from compas_thrust.diagrams.form import _form
from compas_thrust.diagrams.form import remove_feet

from compas_thrust.plotters.plotters import plot_form

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    for j in [1,2]:
        i = 5
        print('\n')
        print('Beginning of Optimisation')

        # file = '/Users/mricardo/compas_dev/me/minmax/radial/0'+ str(j)+ '_0'+ str(i)+ '_sym.json'
        # file_save = '/Users/mricardo/compas_dev/me/minmax/radial/0'+ str(j)+ '_0'+ str(i)+ '_calc.json'
        # file_complete = '/Users/mricardo/compas_dev/me/minmax/radial/0'+ str(j)+ '_0'+ str(i)+ '_complete.json'
        # file_complete_save = '/Users/mricardo/compas_dev/me/minmax/radial/0'+ str(j)+ '_0'+ str(i)+ '_complete_mint.json'

        file = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_TNA.json'
        file_ind = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_ind.json'
        file_save = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_calc.json'
        form = FormDiagram.from_json(file)
        form.plot()
        plot_form(form).show()
        # form_complete = replicate2(form, '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_comp.json', plot=True)
        # form_complete.to_json(file_save)
        form = remove_feet(form, plot = True, openings= 12.0)
        form.plot()
        form.to_json(file_ind)
        form = _form(form, keep_q=True)
        plot_form(form).show()
        form = create_sym_2(form)
        plot_form(form).show()

        # Initial parameters

        tmax = None # form.attributes['tmax']
        bounds_width = 0.1
        use_bounds = True
        qmax = 8.0
        indset = None
        nsol = 9
        sol = 0
        sols = []
        forms = []

        # plot_form(form,radius=0.04).show()

        # Optimisation Routine - Many trials and Many solutions

        for k in range(100):
            print('Optimisation trial {0}.'.format(k))
            if sol == 0 and k % 3 == 0 and k > 0:
                form = _form(form)
                print('Shuffle the form')
                # tmax = 10
            if k > 0 and k % 3 is not 0:
                q = [[form.get_edge_attribute((u,v), 'q')] for u, v in form.edges_where({'is_ind' : True})]
                qi_max  = max(array(q))
                if qi_max > 0.95 * qmax:
                    qmax = qmax * 1.25
                    print('upgrate on qmax')
                else:
                    bounds_width = bounds_width/k
                    use_bounds = True
            fopt, qopt = optimise_single(form, qmax=qmax, solver='devo',
                                            polish='slsqp',
                                            population=800,
                                            generations=600,
                                            printout=50,
                                            tol=0.01,
                                            t = tmax,
                                            opt_max=False,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='loadpath',
                                            indset=indset)
            q = [attr['q'] for u, v, attr in form.edges(True)]
            qmin  = min(array(q))
            # plot_form(form).show()
            if qmin > -0.1: # and check_constraints(form) < 1.0:
                forms.append(deepcopy(form))
                sols.append(fopt)
                form.to_json(file_save)
                print('Optimisation found a result after {0} trials.'.format(k))
                sol = sol + 1
                if sol >= nsol: # or abs((f0 - fopt)/f0) < 0.001:
                    break

        # Printing Results

        print('Solutions for example: {0}'.format(i))
        print(sols)
        n = argmin(sols)
        print('Optimum Solution: {0}'.format(n))
        form = forms[n]
        reactions(form)
        print('Horizontal checks: {0}'.format(horizontal_check(form)))
        check_constraints(form)
        oveview_forces(form)
        form.to_json(file_save)
