from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
# from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.plotters.plotters import plot_form
from compas_thrust.plotters.plotters import plot_force
from compas_thrust.algorithms import initialize_problem
from compas_thrust.algorithms import update_tna

from compas.geometry import is_point_on_segment
from compas.geometry import intersection_segment_segment
from compas.geometry import is_intersection_line_line
from compas.geometry import intersection_line_line
from numpy import shape


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # for i in range(2,9):
    # print('Form: ',str(i))
    # filesym = '/Users/mricardo/compas_dev/me/minmax/radial/02_0'+ str(i) +'_calc.json'
    # form = FormDiagram.from_json(filesym)
    # plot_form(form).show()
    # file = '/Users/mricardo/compas_dev/me/minmax/radial/02_0'+ str(i) +'_complete.json'
    # form = FormDiagram.from_json(file)
    # plot_form(form,radius = 0.05, show_q= False, fix_width=True, max_width=5).show()
    # form, force = update_tna(form)

    # plot_force(force, form).show()
    # plot_form(form, max_width = 5, fix_width=False, show_q= False).show()
    # reactions(form, plot=True)
    # form.to_json(file)

    # lines = [
    #         [   [0.0,0.0,0.0],[0.0,2.0,0.0]     ],
    #         [   [0.0,0.0,0.0],[0.0,-2.0,0.0]    ],
    #         [   [0.0,0.0,0.0],[2.0,0.0,0.0]     ],
    #         [   [0.0,0.0,0.0],[-2.0,0.0,0.0]    ]
    #         ]
    # form = FormDiagram.from_lines(lines, delete_boundary_face=False)
    # form.update_default_vertex_attributes({'is_roller': False})
    # form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    # print(form.number_of_edges())
    # for key in form.vertices():
    #     if form.vertex_coordinates(key) == [0,0,0]:
    #         pass
    #     else:
    #         form.set_vertex_attribute(key, 'is_fixed', True)
    #         form.set_vertex_attribute(key, 'is_anchored', True)
    # form.plot()
    # plot_form(form, fix_width= 0.1, simple=True).show()
    # force = ForceDiagram.from_formdiagram(form)
    # plot_force(force,form)

    file = '/Users/mricardo/compas_dev/me/minmax/radial/02_02_complete.json'
    form = ForceDiagram.from_json(file)
    print(form.number_of_edges())
    # form.plot()
    # plot_form(form, fix_width=True, max_width=5, show_q=False, simple=True).show()
    # plot_form(form).show()
    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = initialize_problem(form)

    print(ind)
    C = C.todense()

    for i in range(C.shape[0]):
        print(C[i,])


