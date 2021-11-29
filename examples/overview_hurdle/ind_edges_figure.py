from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import reciprocal_from_form
from compas_plotters import MeshPlotter
from compas_tno.algorithms import apply_sag
from compas_tno.plotters import plot_independents
from compas_tno.plotters import plot_force_independents

from compas_tno.problems import initialise_form
from compas_tno.algorithms import q_from_qid
from compas_tno.algorithms import xyz_from_q

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


data = {
    'type': 'cross_fd',  # 'cross_with_diagonal' 'cross_fd', 'ortho_fd'
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 4,
    'fix': 'corners',
}

# for bf in [1.0, 2.0, 5.0, 10.0, 20.0]:

for bf in [5.0]:

    form = FormDiagram.from_library(data)

    for key in form.vertices():
        form.vertex_attribute(key, 'pz', -1.0)

    print('Number of edges:', form.number_of_edges())

    apply_sag(form, boundary_force=bf)

    # plotter = MeshPlotter(form, figsize=(8, 8))
    # plotter.draw_edges(width=1.0)
    # plotter.draw_vertices(keys=form.fixed(), facecolor='000000')
    # plotter.show()

    force = reciprocal_from_form(form)

    d = 27.0  # used 27.0 for plot @ overview paper
    number_ind = False

    lines = find_lines(force, d)

    # plotter = MeshPlotter(force, figsize=(8, 8))
    # plotter.draw_edges()
    # plotter.show()

    M = initialise_form(form)

    plot_independents(form, width=1.0, radius=None, number_ind=number_ind).show()
    force_plot = plot_force_independents(force, form, number_ind=number_ind, width=1.0, radius=None)
    force_plot.draw_lines(lines)
    force_plot.show()
    plot_form(form, show_q=True, thick='q').show()

    cut_pattern = []
    ind = []
    i = 0
    q = []
    for edge in form.edges_where({'_is_edge': True}):
        q.append(form.edge_attribute(edge, 'q'))
        if form.edge_attribute(edge, 'is_ind'):
            cut_pattern.append(True)
            ind.append(i)
        i += 1

    q = array(q).reshape(-1, 1)
    q0 = q.copy()

    print(q0[ind])

    for j in range(len(ind)):

        q[ind] = q0[ind]
        q[ind[j]] = q[ind[j]]*1.5
        # q[ind[j]] = q[ind[j]] - 1.0
        print('Modification j=', j)
        print(q[ind])

        q = q_from_qid(q, M.ind, M.Edinv, M.Ei, M.ph)
        M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

        i = 0
        for edge in form.edges_where({'_is_edge': True}):
            form.edge_attribute(edge, 'q', float(q[i]))
            i += 1

        i = 0
        for key in form.vertices():
            form.vertex_attribute(key, 'x', M.X[i, 0])
            form.vertex_attribute(key, 'y', M.X[i, 1])
            form.vertex_attribute(key, 'z', M.X[i, 2])
            i += 1

        force = reciprocal_from_form(form)

        lines = find_lines(force, d)

        plot_independents(form, width=1.0, radius=None).show()
        force_plot = plot_force_independents(force, form, width=1.0, number_ind=number_ind, radius=None)
        force_plot.draw_lines(lines)
        force_plot.show()
        plot_form(form, show_q=True, thick='q').show()

        q[ind[j]] = q[ind[j]]/1.5
        # q[ind[j]] = q[ind[j]] + 1.0

q = q0
print('Plot_normal')
print(q[ind])

q = q_from_qid(q, M.ind, M.Edinv, M.Ei, M.ph)
M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

i = 0
for edge in form.edges_where({'_is_edge': True}):
    form.edge_attribute(edge, 'q', float(q[i]))
    i += 1

i = 0
for key in form.vertices():
    form.vertex_attribute(key, 'x', M.X[i, 0])
    form.vertex_attribute(key, 'y', M.X[i, 1])
    form.vertex_attribute(key, 'z', M.X[i, 2])
    i += 1

force = reciprocal_from_form(form)

lines = find_lines(force, d)

plot_independents(form, width=1.0, radius=None).show()
force_plot = plot_force_independents(force, form, number_ind=number_ind, width=1.0, radius=None)
force_plot.draw_lines(lines)
force_plot.show()
