from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import move_pattern_to_origin

from compas.datastructures import Mesh
from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import normalize_vector

from compas.geometry import Scale
from compas.colors import Color

from numpy import array

import os

def organise_indices(form: FormDiagram):

    edges = []
    for edge in form.edges():
        edges.append(form.edge_coordinates(*edge))

    edges.sort()

    mesh = Mesh.from_lines(edges, delete_boundary_face=True)
    form = FormDiagram.from_mesh(mesh)
    form.set_boundary_supports()
    form.delete_boundary_edges()

    return form

xspan = yspan = [0, 10.0]
thk = 0.50
discretisation = [16, 20]
xc = yc = (xspan[1] - xspan[0])/2
radius = (xspan[1] + xspan[0])/2

dome = Shape.create_dome(thk=thk, radius=radius, t=1.0)

sols = {}

## ni = 2 | ns = 5 in paper MRC
for ni in [4]:
    sols[ni] = {}
    for ns in [2]:
        sols[ni][ns] = {'SWT': None, 'FOBJ': None, 'NIND': None}

        print('\n\n*** ni = {} | ns = {} ***'.format(ni, ns))

        title = 'circular_sigularity_ni_{}_ns_{}.json'.format(ni, ns)
        formpath = '/Users/mricardo/compas_dev/me/pattern/dome/split_support/parametric/form-'+title+'.json'

        form = FormDiagram.from_json(formpath)
        form.parameters['type'] = 'form_sigularity_ni_{}_ns_{}.json'.format(ni, ns)

        uv_i = form.uv_index()
        indices = {edge: uv_i[edge] for edge in form.edges_where({'_is_edge': True})}

        plotter = TNOPlotter(form)
        plotter.draw_form(scale_width=False)
        plotter.draw_edgelabels(text=indices)
        plotter.show()

        # form = organise_indices(form)

        # uv_i = form.uv_index()
        # indices = {edge: uv_i[edge] for edge in form.edges_where({'_is_edge': True})}

        # plotter = TNOPlotter(form)
        # plotter.draw_form(scale_width=False)
        # plotter.draw_edgelabels(text=indices)
        # plotter.show()

        move_pattern_to_origin(form)

        # form.edges_attribute('is_ind', False)
        form.update_indset()

        vector_supports = []
        vectors_plot = []
        base_plot = []

        for key in form.vertices_where({'is_fixed': True}):
            x, y, z = form.vertex_coordinates(key)
            dXbi = [0, 0, 0]

            # dXbi = normalize_vector([(x - xc), (y - yc), 0.0])  # if intended to analyse with outward displacement
            # vectors_plot.append(Vector(*dXbi))
            # base_plot.append(Point(x, y, z))

            if x - xc > 0.1:
                dXbi = [1, 0, 0]
                vectors_plot.append(Vector(*dXbi))
                base_plot.append(Point(x, y, z))
            if x - xc < -0.1:
                dXbi = [-1, 0, 0]
                vectors_plot.append(Vector(*dXbi))
                base_plot.append(Point(x, y, z))

            vector_supports.append(dXbi)

        dXb = array(vector_supports)

        constraints = ['funicular', 'envelope', 'reac_bounds']

        analysis = Analysis.create_compl_energy_analysis(form,
                                                        dome,
                                                        printout=True,
                                                        support_displacement=dXb,
                                                        max_iter=1000,
                                                        solver='IPOPT',
                                                        starting_point='loadpath')

        analysis.optimiser.set_constraints(constraints)
        analysis.optimiser.settings['tol_inds'] = 1e-6
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()

        from compas_tno.problems.problems import plot_svds
        M = analysis.optimiser.M
        plot_svds(M)

        swt = form.lumped_swt()

        print('# inds', form.number_of_independents())

        # pzs = {key: round(form.vertex_attribute(key, 'pz'), 1) for key in form.vertices()}

        # plotter = TNOPlotter(form)
        # plotter.settings['color.edges.independent'] = Color.blue()
        # plotter.draw_form(scale_width=False)
        # plotter.draw_vertexlabels(text= pzs)
        # plotter.show()

        plotter = TNOPlotter(form)
        plotter.settings['color.edges.independent'] = Color.blue()
        plotter.draw_form_independents()
        plotter.draw_supports()
        plotter.show()

        view = Viewer(form, shape=dome)
        view.settings['camera.show.grid'] = False
        view.settings['camera.distance'] = 35
        view.settings['camera.target'] = [5, 5, 0]
        view.settings['camera.rz'] = 45
        view.settings['camera.rx'] = 60
        view.draw_form()
        # view.draw_shape()
        view.show()

        analysis.run()

        folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/dome/split/parametric/', form.parameters['type']) + '/'
        # folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/dome/outwards/parametric/', form.parameters['type']) + '/'
        os.makedirs(folder, exist_ok=True)
        address = folder + title + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
        if analysis.optimiser.exitflag == 0:
            print('Saving form to:', address)
            form.to_json(address)

            sols[ni][ns]['SWT'] = swt
            sols[ni][ns]['FOBJ'] = analysis.optimiser.fopt
            sols[ni][ns]['NIND'] = form.number_of_independents()

plotter = TNOPlotter(form)
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    plotter.draw_vector(vector=vector, base=base)
plotter.draw_form()
plotter.draw_supports()
plotter.draw_cracks()
plotter.show()

view = Viewer(form, shape=dome)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form()
view.draw_cracks()
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()

M = analysis.optimiser.M
from compas_tno.algorithms import equilibrium_residual
equilibrium_residual(M.q, M)

for ni in sols:
    for ns in sols[ni]:
        print(ni, ns, sols[ni][ns]['SWT'], sols[ni][ns]['FOBJ'], sols[ni][ns]['NIND'])
