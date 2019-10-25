
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single
# from compas_thrust.algorithms.mult_inds import optimise_single

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check
from compas_thrust.algorithms.scale import scale_form

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.utilities.constraints import set_height_constraint
from compas_thrust.utilities.constraints import set_cross_vault_heights
from compas_thrust.utilities.constraints import set_dome_heights
from compas_thrust.utilities.constraints import distance_target
from compas_thrust.utilities.constraints import set_cross_vault_heights
from compas_thrust.utilities.constraints import set_pavillion_vault_heights
from compas_thrust.utilities.constraints import circular_heights

from numpy.random import rand
from compas.utilities import geometric_key

from compas_thrust.utilities.loads import set_dome_loads

from compas_thrust.diagrams.form import overview_forces
from compas_thrust.utilities.symmetry import create_sym2
from compas_thrust.utilities.symmetry import replicate2
from compas_thrust.utilities.symmetry import create_sym
from compas_thrust.utilities.symmetry import replicate
from compas_thrust.utilities.symmetry import fix_boundaries_complete
from compas_thrust.utilities.symmetry import fix_mid_complete

from compas_thrust.utilities.loads import not_sym_load
from compas_thrust.utilities.loads import fill_load

from compas_thrust.diagrams.form import _form
from compas_thrust.diagrams.form import remove_feet
from compas_thrust.diagrams.form import delete_boundary_edges

from compas_thrust.diagrams.form import adapt_tna
from compas_thrust.diagrams.form import remove_feet

from compas_thrust.plotters.plotters import plot_form
from compas_viewers.meshviewer import MeshViewer
from compas_thrust.algorithms import z_from_form
from compas_thrust.algorithms import z_update

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

    for j in [1]: # j = 1
        if j == 0 or j == 1:
            shapes = [5] # [2,4,5,7,8,10,11]
        else:
            shapes = [5] #range(2,9)
        # ['B1''B2','B3','C1','C2','C3','C4','C5','D1','D2','D3','D4','D5']
        for i in shapes: # j = 2

            print('\n\n------------------ Form ',str(i),'\n')

            start_time = time.time()

            file_complete = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/01_lp.json'
            file_save = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/minthk/2DArch_Le=2,7_t=20.json'
            # file_complete = '/Users/mricardo/compas_dev/me/convex/4bars/diagram.json'
            # file_complete = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0' + str(j) + '_0' + str(i) + '_complete_paper.json'

            form = FormDiagram.from_json(file_complete)
            # form = _form(form, keep_q=False)
            # form.set_edges_attribute('q', value=3.0)
            # form.set_vertices_attribute('ub', value = 1.0)
            # form.set_vertices_attribute('lb', value = None)
            # plot_form(form, show_q=False, simple=True, max_width=3.0).show()   
            # form = set_pavillion_vault_heights(form, ub_lb=False, thk=2.0, set_heights=True)
            form = circular_heights(form, thk=0.75)
            check_constraints(form, show= True)
            # for key in form.vertices():
            #     if form.get_vertex_attribute(key, 'is_fixed') == False:
            #         form.set_vertex_attribute(key, 'ub', value = 1.0)
            #         form.set_vertex_attribute(key, 'lb', value = 0.5)
            #         print(form.get_vertex_attribute(key,'ub'),form.get_vertex_attribute(key,'lb'))
            # form = fix_boundaries_complete(form)
            # form = fix_mid_complete(form)
            # form.to_json(file_fixed)
            # print(form.vertices_on_boundary())

            # viewer = MeshViewer()
            # viewer.mesh = form
            # viewer.show()

            # for u,v in form.edges():
            #     # tgt = form.get_vertex_attribute(key,'target')
            #     form.set_edge_attribute((u,v), 'q', value=1.0)
            # plot_form(form, show_q=False, simple=True, max_width=3.0).show()

            # overview_forces(form)

            # viewer = MeshViewer()
            # viewer.mesh = form
            # viewer.show()
            
            # for u,v in form.edges():
            #     form.set_edge_attribute((u,v), 'q', value = 1.0)

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

            tmax = 0.05  #form.attributes['tmax'] # 5.0 # form.attributes['tmax']
            bounds_width = None
            use_bounds = False
            qmax = 500.0
            indset = None
            nsol = 5
            sol = 0
            sols = []
            forms = []
            tol = 0.001

            # print('t',tmax)

            # plot_form(form).show()

            # Optimisation Routine - Many trials and Many solutions

            fopt = 100

            for k in range(5):
                # form = _form(form, keep_q=True)
                fopt, qopt = optimise_single(form, qmax=qmax, solver=None,
                                                polish='slsqp',
                                                printout=100,
                                                tol=0.01,
                                                t = tmax,
                                                opt_max=False,
                                                tension=False,
                                                use_bounds = use_bounds,
                                                bounds_width = bounds_width,
                                                objective='min',
                                                indset=indset,
                                                buttress=True)
                q = array([attr['q'] for u, v, attr in form.edges(True)])
                print('Result iteration {0}: {1}'.format(k,fopt))
                qmin  = min(q)
                if qmin > -0.1 and fopt is not None and check_constraints(form, show= False) < 1.0:
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
                overview_forces(form)

                time_i = end_time - start_time
                ts.append(time_i)
                print('Times elapsed: {0:.3f} sec'.format(time_i))
                print('Number of edges: {0}'.format(form.number_of_edges()))
                ns.append(form.number_of_edges())
                # plot_form(form).show()
                # form.to_json(file_save)

            else:
                print('\n\nFailed on', i)

    print('Final LP, TIMES, N_EDGES, OBJ_FUNC')
    print(lps)
    print(ts)
    print(ns)
    print(obj)

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()

