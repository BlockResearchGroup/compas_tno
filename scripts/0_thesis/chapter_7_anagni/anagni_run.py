from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_selfweight_from_thrust
from compas_tno.problems.problems import initialise_form, plot_svds
from compas.datastructures import Mesh

from compas.colors import Color
import json
import os

thk = 0.25
additional_thk = 0.00
error = 0.0
ro_fill = 14.0

thicknessess = []
t_over_W = []

# Reassuress from Rhinoceros
xspan = [-4.431, 1.313]
yspan = [14.968, 18.312]

for diagram_name in ['mesh-D3']:

    print('Mesh:', diagram_name)

    fill_loads = False
    folder = '/Users/mricardo/compas_dev/me/anagni/'

    x0 = 0.0
    xf = xspan[1] - xspan[0]
    y0 = 0.0
    yf = yspan[1] - yspan[0]
    xc = (xf + x0)/2
    yc = (yf + y0)/2
    discretisation = 12

    corners = [[x0, y0], [xf, y0], [xf, yf], [x0, yf]]
    print('Corners:', corners)

    # ------------------- FormDiagram --------------------

    file_formdiagram = os.path.join(folder, 'meshes', 'CISM', diagram_name + '.json')

    mesh = Mesh.from_json(file_formdiagram)

    vertices, faces = mesh.to_vertices_and_faces()
    form = FormDiagram.from_vertices_and_faces(vertices, faces)

    # path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-D3/sangelo_vault_top_final_mesh-D3_with_fill_n.json'
    # # path = '/Users/mricardo/compas_dev/me/anagni/meshes/CISM/form-D_with_inds.json'
    # form = FormDiagram.from_json(path)

    # form = FormDiagram.create_parametric_form(discretisation=discretisation, lambd=1.0)

    move_pattern_to_origin(form, corners)  # move pattern to origin and find supports on corners

    M = initialise_form(form, printout=True)
    # M = initialise_form(form, printout=True)
    plot_svds(M)


    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False, color=Color.black())
    # plotter.draw_lines(lines)
    plotter.draw_supports(color=Color.red())
    plotter.zoom_extents()  # why zoom extents does not work?
    plotter.show()

    # ----------------------- Point Cloud -----------------------

    file_name = 'sangelo_vault_top_final'  # This point cloud considers the extrados with 25 cm only
    # Note: The fill loads must be added
    pointcloud = folder + file_name + '.json'

    print('Point-cloud_location:', pointcloud)

    points = {'UB': [], 'LB': [], 'FILL': []}

    tol = 0.001

    with open(pointcloud) as json_file:
        data = json.load(json_file)

    for tp in ['UB', 'LB', 'FILL']:
        for key, pt in data[tp].items():
            x, y, z = pt
            pt = [x - xspan[0], y - yspan[0], z]

            if abs(pt[0] - x0) < tol:
                pt[0] = pt[0] - tol * 2
            elif abs(pt[0] - xf) < tol:
                pt[0] = pt[0] + tol * 2

            if abs(pt[1] - y0) < tol:
                pt[1] = pt[1] - tol * 2
            elif abs(pt[1] - yf) < tol:
                pt[1] = pt[1] + tol * 2

            points[tp].append(pt)


    # ------- Create shape given a topology and a point cloud --------

    vault = Shape.from_pointcloud_and_topology(form, points['LB'], points['UB'], fill_pts=points['FILL'], data={'type': 'general', 't': 0.0, 'thk': thk + additional_thk})
    vault.store_normals(plot=False)

    area = vault.middle.area()
    swt = vault.compute_selfweight()

    print('Interpolated Volume Data:')
    print('Self-weight is: {0:.2f}'.format(swt))
    print('Area is: {0:.2f}'.format(area))

    apply_selfweight_from_shape(form, vault)
    apply_envelope_from_shape(form, vault)

    pzt0 = 0
    for key in form.vertices():
        pzi = form.vertex_attribute(key, 'pz')
        pzt0 += pzi

    print('Selfweight NOT considering the fill is:', pzt0)

    # ----- Apply fill loads and bounds

    if fill_loads:
        pzfilldic = {}
        at = 0.0
        for key in vault.fill.vertices():
            pz = form.vertex_attribute(key, 'pz')
            zub = form.vertex_attribute(key, 'ub')
            z_fill = vault.fill.vertex_attribute(key, 'z')
            z_diff = (z_fill - zub)
            if z_diff > 0.01:
                ai = form.vertex_area(key)  # since the form is planar at this point this corresponds to the projected area of the node
                pz_fill = -1 * ai * z_diff * ro_fill  # downwards loads are negative
                pzt = pz + pz_fill
                pzfilldic[key] = round(pz_fill, 1)
                form.vertex_attribute(key, 'pz', pzt)
                # print('improve height of:', key)
                at += ai

    analysis: Analysis = Analysis.create_minthrust_analysis(form, vault,
                                                            printout=True,
                                                            solver='IPOPT',
                                                            max_iter=2000,
                                                            plot=True)
    # analysis.optimiser.set_additional_options(solver_convex='MATLAB')
    analysis.set_up_optimiser()
    analysis.run()

    pzt = form.lumped_swt()
    thrust = form.thrust()
    t_over_wi = thrust/pzt
    print('T/W:', t_over_wi)

    view = Viewer(form, vault)
    view.settings['camera.target'] = [2.8, 1.6, 12]
    view.settings['camera.distance'] = 18
    view.settings['scale.reactions'] = 0.005 * 4
    view.settings['scale.loads'] = 0.005 * 4 * 10
    view.settings['opacity.shapes'] = 0.3
    view.scale_edge_thickness(3.0)
    view.draw_form()
    view.draw_shape()
    view.draw_cracks()
    view.show()

