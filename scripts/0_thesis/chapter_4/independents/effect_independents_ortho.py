from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
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


form = FormDiagram.create_ortho_form(discretisation=6, fix='all')
# form = FormDiagram.create_cross_form(discretisation=6, fix='corners')

apply_selfweight_from_thrust(form)

initialize_tna(form)
# initialize_loadpath(form)

viewer = Viewer(form, show_grid=False)
viewer.settings['scale.edge.thk_absolute'] = 1.5
viewer.settings['scale.reactions'] = 0.02
viewer.settings['camera.target'] = [5, 5, 0]
viewer.draw_thrust(absolute_scale=True)
viewer.draw_thrustsurface()
viewer.draw_reactions()
viewer.show()

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

force = reciprocal_from_form(form)

d = 10.0  # used 27.0 for plot @ overview paper
number_ind = False

M = initialise_form(form)

scale_force = None

plotter = TNOPlotter(form, force=force)
plotter.draw_form_independents()
plotter.draw_supports(color=Color.red())
plotter.draw_force(scale=scale_force, show_independents=True)
plotter.show()

update_geometry(form, M)

viewer = Viewer(form, show_grid=False)
viewer.settings['scale.edge.thk_absolute'] = 1.5
viewer.settings['scale.reactions'] = 0.02
viewer.settings['camera.target'] = [5, 5, 0]
viewer.draw_thrust(absolute_scale=True)
viewer.draw_thrustsurface()
viewer.draw_reactions()
viewer.show()

q0 = M.q.copy()

print(q0[M.ind])

for j in range(len(M.ind)):

    mult = 2.0

    M.q[M.ind] = q0[M.ind]
    M.q[M.ind[j]] = M.q[M.ind[j]]*mult
    # q[ind[j]] = q[ind[j]] - 1.0
    print('Modification j=', j)
    print(M.q[M.ind])

    update_geometry(form, M)

    force = reciprocal_from_form(form)

    save_img = '/Users/mricardo/Documents/ETH/Thesis/PhD_Thesis/figures/chapter_4/source/form_forces/' + 'vault-' + str(j) + '.pdf'
    plotter = TNOPlotter(form, force=force)
    plotter.draw_form_independents()
    plotter.draw_supports()
    plotter.draw_force(scale=scale_force, show_independents=True)
    # plotter.save(save_img)
    plotter.show()

    M.q[M.ind[j]] = M.q[M.ind[j]]/mult

    viewer = Viewer(form, show_grid=False)
    viewer.settings['scale.edge.thk_absolute'] = 1.5
    viewer.settings['scale.reactions'] = 0.02
    viewer.settings['camera.target'] = [5, 5, 0]
    viewer.draw_thrust(absolute_scale=True)
    viewer.draw_thrustsurface()
    viewer.draw_reactions()
    viewer.show()
