
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.utilities.symmetry import create_sym2
from compas_thrust.utilities.symmetry import replicate2
from compas_thrust.utilities.symmetry import create_sym
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import _form
from compas_thrust.diagrams.form import remove_feet
from compas_thrust.diagrams.form import oveview_forces

from compas_thrust.diagrams.form import adapt_tna
from compas_thrust.diagrams.form import adapt_objective
from compas_thrust.diagrams.form import remove_feet

from compas_thrust.plotters.plotters import plot_form

from compas_viewers.meshviewer import MeshViewer

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # for j in [2, 1]:
    #     for i in [5,6,7]:
    print('\n')
    print('Beginning of Optimisation')
    # print('Form: '+str(j)+' Discretization: '+str(i))

    # file = '/Users/mricardo/compas_dev/me/minmax/fan/0'+ str(j)+ '_0'+ str(i)+ '_sym.json'
    # file_save = '/Users/mricardo/compas_dev/me/minmax/fan/0'+ str(j)+ '_0'+ str(i)+ '_calc.json'
    # file_complete = '/Users/mricardo/compas_dev/me/minmax/fan/0'+ str(j)+ '_0'+ str(i)+ '_complete.json'
    # file_complete_save = '/Users/mricardo/compas_dev/me/minmax/fan/0'+ str(j)+ '_0'+ str(i)+ '_complete_min.json'

    file = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_comp.json'
    file_scaled = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_scaled.json'
    file_sym = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_sym.json'
    file_save = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_calc.json'
    file_complete_save = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_comp_calc.json'

    # form = FormDiagram.from_json(file_scaled)
    # form = _form(form)
    # form.plot()
    # plot_form(form).show()

    # adapt_objective(form, objective = 'loadpath', plot = False, delete_face = True)
    # form.to_json(file_scaled)
    
    # form.plot()
    # form = adapt_tna(form, zmax = 7.0, plot = False, delete_face = True)
    # plot_form(form).show()
    # form_complete = replicate2(form, '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_comp.json', plot=True)
    # form_complete.to_json(file_save)
    
    # form = remove_feet(form, plot = False, openings = 12.0)
    
    # form.plot()
    # form = _form(form, keep_q=True)
    # plot_form(form).show()
    # form = create_sym_2(form)
    # plot_form(form).show()

    # form = create_sym(form)
    # form.to_json(file_sym)

    # Initial parameters

    # tmax = None #3.0 # form.attributes['tmax']
    # bounds_width = 2
    # use_bounds = True
    # qmax = 10
    # indset = None
    # nsol = 4
    # sol = 0
    # sols = []
    # forms = []

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
                                        population=300,
                                        generations=200,
                                        printout=None,
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
        if qmin > -0.1 and check_constraints(form, show= False) < 0.1:
            forms.append(deepcopy(form))
            sols.append(fopt)
            form.to_json(file_save)
            print('Optimisation found a result after {0} trials.'.format(k))
            sol = sol + 1
            if sol >= nsol: # or abs((f0 - fopt)/f0) < 0.001:
                break

    # Printing Results

    print('Solutions for example:') # {0}'.format(i))
    print(sols)
    n = argmin(sols)
    print('Optimum Solution: {0}'.format(n))
    form = forms[n]
    reactions(form)
    print('Horizontal checks: {0}'.format(horizontal_check(form)))
    check_constraints(form)
    oveview_forces(form)
    form.to_json(file_save)

    # Replicate and Save Complete

    form = FormDiagram.from_json(file_save)
    oveview_forces(form)
    plot_form(form).show()
    form_ = replicate(form, file, plot=True)
    reactions(form_)
    check_constraints(form_)
    form_.to_json(file_complete_save)

