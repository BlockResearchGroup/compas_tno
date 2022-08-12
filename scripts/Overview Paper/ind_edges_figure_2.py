from compas_tno import viewers
from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import reciprocal_from_form
from compas_plotters import MeshPlotter
from compas_tno.algorithms import apply_sag
from compas_tno.plotters import plot_independents
from compas_tno.plotters import plot_force_independents
import compas_tno

from compas_tno.problems import initialise_form
from compas_tno.algorithms import q_from_qid
from compas_tno.algorithms import xyz_from_q

from compas_tno.viewers import Viewer

from compas.datastructures import mesh_weld

from numpy import array


def find_lines(force, d):
    [xc, yc, _] = force.centroid()

    p1 = [xc - d, yc + d]
    p2 = [xc + d, yc + d]
    p3 = [xc + d, yc - d]
    p4 = [xc - d, yc - d]

    l1 = {'start': p1, 'end': p2}
    l2 = {'start': p2, 'end': p3}
    l3 = {'start': p3, 'end': p4}
    l4 = {'start': p4, 'end': p1}

    return [l1, l2, l3, l4]


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


data = {
    'type': 'ortho_fd',  # 'cross_with_diagonal' 'cross_fd', 'ortho_fd'
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 6,
    'fix': 'all',
}

# for bf in [1.0, 2.0, 5.0, 10.0, 20.0]:

for bf in [5.0]:

    form = FormDiagram.from_library(data)
    # form = mesh_weld(form)

    for edge in form.edges_on_boundary():
        form.edge_attribute(edge, '_is_edge', False)

    for key in form.vertices():
        form.vertex_attribute(key, 'pz', 1.0)

    ad = compas_tno.get('test.json')
    form.to_json(ad)
    print(ad)

    print('Number of edges:', form.number_of_edges())

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

    lines = find_lines(force, d)

    # plotter = Plotter(force, figsize=(8, 8))
    # plotter.draw_edges()
    # plotter.show()

    M = initialise_form(form)

    plot_independents(form, width=1.0, radius=None, number_ind=number_ind).show()
    force_plot = plot_force_independents(force, form, number_ind=number_ind, width=1.0, radius=None)
    force_plot.draw_lines(lines)
    force_plot.show()
    # plot_form(form, show_q=True, thick='q').show()

    update_geometry(form, M)

    address = '/Users/mricardo/compas_dev/compas_tno/data/form_q={}.json'.format(0)
    form.to_json(address)

    viewer = Viewer(form)
    viewer.draw_thrust()
    viewer.show()

    q0 = M.q.copy()

    print(q0[M.ind])

    for j in [2, 6]:  # range(len(M.ind)):

        mult = 4.0

        M.q[M.ind] = q0[M.ind]
        M.q[M.ind[j]] = M.q[M.ind[j]]*mult
        # q[ind[j]] = q[ind[j]] - 1.0
        print('Modification j=', j)
        print(M.q[M.ind])

        update_geometry(form, M)

        address = '/Users/mricardo/compas_dev/compas_tno/data/form_q={}.json'.format(j)
        form.to_json(address)

        force = reciprocal_from_form(form)

        lines = find_lines(force, d)

        # plot_independents(form, width=1.0, radius=None).show()
        force_plot = plot_force_independents(force, form, width=1.0, number_ind=number_ind, radius=None)
        force_plot.draw_lines(lines)
        force_plot.show()
        # plot_form(form, show_q=True, thick='q').show()

        M.q[M.ind[j]] = M.q[M.ind[j]]/mult
        # M.q[ind[j]] = M.q[ind[j]] + 1.0

        viewer = Viewer(form)
        viewer.draw_thrust()
        viewer.show()


M.q = q0
print('Plot_normal')
print(M.q[M.ind])

update_geometry(form, M)

force = reciprocal_from_form(form)

lines = find_lines(force, d)

plot_independents(form, width=1.0, radius=None).show()
force_plot = plot_force_independents(force, form, number_ind=number_ind, width=1.0, radius=None)
force_plot.draw_lines(lines)
force_plot.show()

viewer = Viewer(form)
viewer.settings['size.edge.max_thickness'] = viewer.settings['size.edge.max_thickness']/4.0
viewer.draw_thrust()
viewer.show()
