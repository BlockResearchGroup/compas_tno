
from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.ind_based import optimise_single
# from compas_tno.algorithms.mult_inds import optimise_single

from compas_tno.algorithms.equilibrium import reactions
from compas_tno.algorithms.equilibrium import horizontal_check
from compas_tno.algorithms.scale import scale_form

from compas_tno.utilities.constraints import check_constraints
from compas_tno.utilities.constraints import set_height_constraint
from compas_tno.utilities.constraints import set_cross_vault_heights
from compas_tno.utilities.constraints import set_dome_heights
from compas_tno.utilities.constraints import distance_target
from compas_tno.utilities.constraints import set_cross_vault_heights
from compas_tno.utilities.constraints import set_pavillion_vault_heights

from numpy.random import rand
from compas.utilities import geometric_key

from compas_tno.utilities.loads import set_dome_loads

from compas_tno.diagrams.form import overview_forces
from compas_tno.utilities.symmetry import create_sym2
from compas_tno.utilities.symmetry import replicate2
from compas_tno.utilities.symmetry import create_sym
from compas_tno.utilities.symmetry import replicate
from compas_tno.utilities.symmetry import fix_boundaries_complete
from compas_tno.utilities.symmetry import fix_mid_complete

from compas_tno.utilities.loads import not_sym_load
from compas_tno.utilities.loads import fill_load

from compas_tno.diagrams.form import _form
from compas_tno.diagrams.form import remove_feet
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.diagrams.form import adapt_tna
from compas_tno.diagrams.form import remove_feet

from compas_tno.plotters import plot_form
from compas_viewers.meshviewer import MeshViewer
from compas_tno.algorithms import z_from_form
from compas_tno.algorithms import z_update

from copy import deepcopy
from numpy import array
from numpy import argmin

import time

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    lps = []
    ts = []
    ns = []
    obj = []

    for j in [2]: # j = 1
        if j == 0 or j == 1:
            shapes = [11] #[2,4,5,7,8,10,11]
        else:
            shapes = [5] #range(2,9)
        # ['B1''B2','B3','C1','C2','C3','C4','C5','D1','D2','D3','D4','D5']
        for i in shapes: # j = 2

            print('\n\n------------------ Form ',str(i),'\n')

            start_time = time.time()

            # Files for calculating the square very fast
            # file = '/Users/mricardo/compas_dev/me/loadpath/Corner/discretize/0'+str(j)+'_0'+str(i)+'_complete.json'
            # file_save = '/Users/mricardo/compas_dev/me/loadpath/Corner/discretize/0'+str(j)+'_0'+str(i)+'_fit_sixpartite.json'
            file_complete = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0' + str(j) + '_0' + str(i) + '_complete_paper.json'

            # file_save = '/Users/mricardo/compas_dev/me/loadpath/Corner/topology/'+i+'_lp.json'
            # file_complete = '/Users/mricardo/compas_dev/me/loadpath/Corner/topology/'+i+'_complete.json'

            form = FormDiagram.from_json(file_complete)
            # plot_form(form, show_q=False, simple=True, max_width=3.0).show()
            # form = set_pavillion_vault_heights(form, ub_lb=False, thk=2.0, set_heights=True)
            # form = fix_boundaries_complete(form)
            # form = fix_mid_complete(form)
            # form.to_json(file_fixed)
            plot_form(form, show_q=False, max_width=2, simple=True, radius=0.02).show()
            # print(form.vertices_on_boundary())

            # viewer = MeshViewer()
            # viewer.mesh = form
            # viewer.show()

            # for u,v in form.edges():
            #     # tgt = form.vertex_attribute(key,'target')
            #     form.edge_attribute((u,v), 'q', value=1.0)
            # plot_form(form, show_q=False, simple=True, max_width=3.0).show()

            # overview_forces(form)

            # viewer = MeshViewer()
            # viewer.mesh = form
            # viewer.show()

            # for u,v in form.edges():
            #     form.edge_attribute((u,v), 'q', value = 1.0)

            # form = create_sym(form, keep_q= True)
            # form = delete_boundary_edges(form)
            # form = set_cross_vault_heights(form, weights=True)
            # form = set_pavillion_vault_heights(form)
            # form = set_dome_heights(form)
            # plot_form(form, show_q=True).show()
            # form.to_json(file)

            # plot_form(form, max_width=2.0, show_q=True).show()

            # overview_forces(form)
            # form = set_height_constraint(form)
            # z = [form.vertex_coordinates(key)[2] for key in form.vertices()]
            # zmax = max(z)
            # r = zmax/1.0
            # form = scale_form(form, r)
            # form = not_sym_load(form)
            # form = fill_load(form)
            # plot_form(form).show()
            # plot_form(form, fix_width=True, max_width=1.0).show()
            # viewer = MeshViewer()
            # viewer.mesh = form
            # viewer.show()
            # overview_forces(form)
            # print('Horizontal checks: {0}'.format(horizontal_check(form)))
            # check_constraints(form, show=True)

            # form.to_json(file)

            # Initial parameters

            tmax = None # 5.0 # form.attributes['tmax']
            bounds_width = 5.0
            use_bounds = False
            qmax = 100.0
            indset = None
            nsol = 5
            sol = 0
            sols = []
            forms = []
            tol = 0.001
            ni = []
            # indset = '0.000,1.563,0.000,2.813,2.461,0.000,8.437,1.563,0.000,5.000,2.187,0.000,10.000,7.813,0.000,10.000,1.563,0.000,7.812,1.367,0.000,7.812,7.812,0.000,8.945,2.812,0.000,1.914,7.813,0.000,2.188,7.812,0.000,8.242,7.188,0.000,1.055,2.812,0.000,7.812,0.820,0.000,0.547,2.188,0.000,5.312,0.000,0.000,4.688,0.000,0.000,3.438,7.852,0.000,8.984,4.063,0.000,4.063,6.445,0.000,9.062,9.648,0.000,0.938,9.297,0.000,0.703,0.938,0.000,5.937,5.000,0.000,2.539,5.937,0.000,0.430,6.562,0.000,5.312,9.414,0.000,5.000,9.063,0.000,8.437,0.781,0.000,5.938,10.000,0.000,8.594,2.813,0.000,0.000,5.312,0.000,7.813,0.273,0.000,4.688,8.828,0.000,9.180,0.938,0.000,6.484,5.312,0.000,1.406,7.187,0.000,4.062,3.047,0.000,9.727,7.813,0.000,6.562,8.281,0.000,8.437,0.391,0.000,2.813,10.000,0.000'
            # q = array([attr['q'] for u, v, attr in form.edges(True)])
            # ind = []
            # for u, v in form.edges():
            #     if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
            #         ind.append(uv_i[(u, v)])
            # print(ind)
            # test = rand(len(ind)) * qmax
            # q[ind] = test

            # print( hdsiodh)

            # plot_form(form,radius=0.04).show()

            # Optimisation Routine - Many trials and Many solutions

            fopt = 100

            for k in range(5):
                # form = fill_load(form)
                form = _form(form, keep_q=True)
                print('Optimisation trial {0}.'.format(k))
                if sol == 0 and k % 4 == 0 and k > 0:
                    form = _form(form, keep_q=False)
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
                                                printout=100,
                                                tol=0.01,
                                                t = tmax,
                                                opt_max=False,
                                                tension=False,
                                                use_bounds = use_bounds,
                                                bounds_width = bounds_width,
                                                objective='loadpath',
                                                indset=indset,
                                                buttress=False)
                # plot_form(form, show_q=False).show()
                q = array([attr['q'] for u, v, attr in form.edges(True)])
                print('Result iteration {0}: {1}'.format(k,fopt))
                qmin  = min(q)
                if qmin > -0.1 and fopt is not None and fopt < 1000.0 and check_constraints(form, show= False) < 1.0:
                    if sols == []:
                        end_time = time.time() # take time of first solution only
                    forms.append(deepcopy(form))
                    sols.append(fopt)
                    print('Optimisation found a result after {0} trials.'.format(k+1))
                    sol = sol + 1
                    if (sol >= nsol or abs(fopt-sols[len(sols)-2]) < 0.001) and len(sols) > 1:
                        break

            # Printing Results

            if sols is not []:

                print('Solutions for example: ')
                print(sols)
                n = argmin(sols)
                print('Optimum Solution: {0}'.format(n))
                form = forms[n]
                lps.append(form.attributes['loadpath'])
                obj.append(sols[n])
                # reactions(form)
                # print('Horizontal checks: {0}'.format(horizontal_check(form)))
                # overview_forces(form)

                time_i = end_time - start_time
                ts.append(time_i)
                print('Times elapsed: {0:.3f} sec'.format(time_i))
                print('Number of edges: {0}'.format(form.number_of_edges()))
                ns.append(form.number_of_edges())
                plot_form(form).show()
                form.to_json(file_save)

            else:
                print('\n\nFailed on', i)

            # form_old = FormDiagram.from_json(file_save)
            # if form.attributes['loadpath'] < form_old.attributes['loadpath']:
            #     form.to_json(file_save)
            # else:
            #     print('The optimisation before: {0} was better than now: {1}'.format(form_old.attributes['loadpath'],form.attributes['loadpath']))
            # print(lps)
            # print(obj)
            # print(ts)
            # print(ns)

            # Replicate and Save Complete

            # form = FormDiagram.from_json(file_complete)
            # form_ = replicate(form, file_complete, plot=True)
            # reactions(form_)
            # overview_forces(form_)
            # check_constraints(form_)
            # form_.to_json(file_save)
            # form = form_

    print('Final LP, TIMES, N_EDGES, OBJ_FUNC')
    print(lps)
    print(ts)
    print(ns)
    print(obj)

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()

