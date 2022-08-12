from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis

from compas_tno.algorithms import reciprocal_from_form
from compas_tno.algorithms import apply_sag
from compas_tno.algorithms import q_from_qid
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import initialise_form
from compas_tno.problems import initialize_tna
from compas_tno.plotters import TNOPlotter
from compas_tno.problems.initialize import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_thrust
from compas_tno.viewers import Viewer

from compas_plotters import Plotter
from compas.geometry import Line
from compas.geometry import Point
from compas.colors import Color

from numpy import array

def update_geometry(form, M):

    M.q = q_from_qid(M.q, M.ind, M.Edinv, M.Ei, M.ph)
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    i = 0
    for edge in form.edges_where({'_is_edge': True}):
        form.edge_attribute(edge, 'q', float(M.q[i]))
        i += 1

    i = 0
    for key in form.vertices():
        form.vertex_attribute(key, 'x', M.X[i, 0])
        form.vertex_attribute(key, 'y', M.X[i, 1])
        form.vertex_attribute(key, 'z', M.X[i, 2])
        i += 1


form = FormDiagram.create_cross_form(discretisation=6, fix='corners')
shape = Shape.create_crossvault()

analysis: Analysis = Analysis.create_minthk_analysis(form, shape)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

M = analysis.optimiser.M

# initialize_tna(form)
# initialize_loadpath(form)

# for edge in form.edges_where({'_is_edge': True}):
#     q = form.edge_attribute(edge, 'q')
#     form.edge_attribute(edge, 'q', 2*q)

# initialize_fdm(form)

print('Number of edges:', form.number_of_real_edges())

# apply_sag(form, boundary_force=bf)

# plotter = MeshPlotter(form, figsize=(8, 8))
# plotter.draw_edges(text={edge: str(form.edge_attribute(edge, '_is_edge')) for edge in form.edges()}, width=1.0)
# plotter.draw_vertices(keys=form.fixed(), facecolor='FF0000')
# print(list(form.faces()))
# plotter.draw_faces(text={fkey: str(fkey) for fkey in form.faces()})
# plotter.show()

force = reciprocal_from_form(form, plot=False, restore_form_topology=True)

d = 10.0  # used 27.0 for plot @ overview paper
number_ind = False

scale_force = 0.030

# plotter = TNOPlotter(form, force=force)
# # plotter.draw_form_independents()
# # plotter.draw_form(scale_width=False)
# # plotter.draw_supports(color=Color.red())
# # plotter.draw_force(scale=scale_force)
# plotter.draw_force(scale=scale_force, show_independents=True)
# plotter.app.add(Point(25.0, 5.0))
# plotter.app.add(Point(25.0, 13.0))
# plotter.app.add(Point(25.0, -3.0))
# plotter.app.add(Point(9.0, 5.0))
# plotter.app.add(Point(9.0, 13.0))
# plotter.app.add(Point(9.0, -3.0))
# plotter.show()

update_geometry(form, M)

viewer = Viewer(form, show_grid=False)
viewer.settings['scale.edge.thk_absolute'] = 1.25
viewer.settings['scale.reactions'] = 0.004  # divided by 10 wit regards to the continuously suupported
viewer.settings['camera.target'] = [5, 5, 0]
viewer.draw_thrust(absolute_scale=True)
viewer.draw_thrustsurface()
viewer.draw_reactions()
viewer.show()

q0 = M.q.copy()

print(q0[M.ind])

for j in range(len(M.ind)):

    if j in [0, 2, 6, 7]:
        mult = 2.5
    else:
        mult = 1.6

    M.q[M.ind] = q0[M.ind]
    M.q[M.ind[j]] = M.q[M.ind[j]]*mult
    # q[ind[j]] = q[ind[j]] - 1.0
    print('Modification j=', j)
    print(M.q[M.ind])

    update_geometry(form, M)

    force = reciprocal_from_form(form)

    # save_img = '/Users/mricardo/Documents/ETH/Thesis/PhD_Thesis/figures/chapter_4/source/form_forces/' + 'vault-' + str(j) + '.pdf'
    # plotter = TNOPlotter(form, force=force)
    # # plotter.draw_form_independents()
    # # plotter.draw_form(scale_width=False)
    # # plotter.draw_supports(color=Color.red())
    # # plotter.draw_force(scale=scale_force)
    # plotter.draw_force(scale=scale_force, show_independents=True)
    # plotter.app.add(Point(25.0, 5.0))
    # plotter.app.add(Point(25.0, 13.0))
    # plotter.app.add(Point(25.0, -3.0))
    # plotter.app.add(Point(9.0, 5.0))
    # plotter.app.add(Point(9.0, 13.0))
    # plotter.app.add(Point(9.0, -3.0))
    # plotter.show()

    M.q[M.ind[j]] = M.q[M.ind[j]]/mult

    viewer = Viewer(form, show_grid=False)
    viewer.settings['scale.edge.thk_absolute'] = 1.25
    viewer.settings['scale.reactions'] = 0.004  # divided by 10 wit regards to the continuously suupported
    viewer.settings['camera.target'] = [5, 5, 0]
    viewer.draw_thrust(absolute_scale=True)
    viewer.draw_thrustsurface()
    viewer.draw_reactions()
    viewer.show()
