from compas_plotters import MeshPlotter
import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_mesh
from compas_tno.viewers import view_meshes
from compas_tno.viewers import view_solution
from compas_plotters import MeshPlotter
import math

from numpy import linspace
from compas.utilities import geometric_key
from compas_tno.datastructures import MeshDos

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN THICKNESS FOR POINTED ARCH --------------------
# ----------------------------------------------------------------------


sols = {}
loads = {}

Rs = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]
degs = [20]

# FAN
# degs = [0]
# Rs = [11.4027]

# degs = [10]
# Rs = [9.7364]

# degs = [20]
# Rs = [7.9422]

# degs = [30]
# Rs = [6.8896]

# degs = [40]
# Rs = [6.2453]

#CROSS
# degs = [0]
# Rs = [7.2299]

# degs = [10]
# Rs = [6.7900]

# degs = [20]
# Rs = [6.1128257]

# degs = [30]
# Rs = [5.6437]

# degs = [40]
# Rs = [5.3306]

# NON OPTIMUM
# Rs = [5.0, 7.5, 10.0]
# degs = [0, 20, 40]

Rs = [6.5]
degs = [20]

save = True
plot = True
arch_diagonal = False
mesh_diagonal = False


discretisation = 14
thk = 0.50
b = 0.5  # Out of plane dimension  of arch
t = 0.0
#-----------------
type_diagram_loads = 'topology-crossbraced'
type_structure = 'pointed_arch'
type_formdiagram = 'pointed_arch'
type_formdiagram_pattern = 'topology-crosssmooth2'
#-----------------

L = 10.0
x0 = 0.0

# ----------------------- 1. Create Form Diagram (Line) ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'total_nodes': discretisation + 1,
    'x0': x0,
    'L': L,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)

if arch_diagonal:
    total_nodes = discretisation + 1
    x = linspace(x0, x0 + L, num=total_nodes, endpoint=True)  # Continue this
    lines = []
    gkey_fix = []

    for i in range(total_nodes-1):
        xi = x[i]
        yi = xi
        xf = x[i + 1]
        yf = xf
        lines.append([[xi, yi, 0.0], [xf, yf, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, yi, 0.0], precision=6))
        elif i == total_nodes - 2:
            gkey_fix.append(geometric_key([xf, xf, 0.0], precision=6))

    mesh = MeshDos.from_lines(lines)
    form = FormDiagram.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    # plotter = MeshPlotter(form)
    # plotter.draw_edges()
    # plotter.draw_vertices(text={key: key for key in form.vertices()})
    # plotter.show()

for deg in degs:

    sols[deg] = {}

    for R in Rs:

        A = L/(2*R*(math.cos(math.radians(deg)) - 1) + L)
        L_shape = A * L
        x0_shape = x0 + L/2*(1 - A)
        hc_shape = A*math.sqrt(R**2 - (R - L/2)**2)
        print('Shape parameters R A L x0 hc:', R, A, L_shape, x0_shape, hc_shape)

        # ----------------------- 1. Create Arch shape ---------------------------

        data_shape = {
            'type': type_structure,
            'hc': hc_shape,
            'L': L_shape,
            'thk': thk,
            'discretisation': discretisation + 1,
            'b': b,
            't': t,
            'x0': x0_shape
        }

        arch = Shape.from_library(data_shape)
        gradients = True

        # ----------------------- 1. Create fan diagram and 3D shape ---------------------------

        xy_span_shape = [[-L/2*(A - 1), L*(1 + (A - 1)/2)], [-L/2*(A - 1), L*(1 + (A - 1)/2)]]
        data_shape_3D = {
                'type': 'pointed_crossvault',
                'thk': thk,
                'discretisation': discretisation,
                'xy_span': xy_span_shape,
                't': 0.0,
                'hc': hc_shape,
                'hm': None,
                'he': None,
            }
        vault = Shape.from_library(data_shape_3D)
        # view_shapes(vault).show()

        if '-' in [char for char in type_diagram_loads]:
            type_topology = type_diagram_loads.split('-')[-1]
        else:
            type_topology = False

        if type_topology:
            type_diagram_loads = 'topology-' + type_topology
            fd_mesh = 'FormDiagram-' + type_topology
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', 'pointed_crossvault', type_diagram_loads)
            form_3D = FormDiagram.from_json((os.path.join(folder, fd_mesh + '.json')))
        else:
            data_diagram_3D = {
                'type': type_diagram_loads,
                'xy_span': [[0, L], [0, L]],
                'discretisation': discretisation,
                'fix': 'corners',
            }
            form_3D = FormDiagram.from_library(data_diagram_3D)

        if type_formdiagram_pattern:
            file_pattern = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', 'pointed_crossvault', type_formdiagram_pattern, 'FormDiagram-' + type_formdiagram_pattern.split('-')[-1] + '.json')
            form_pattern = FormDiagram.from_json(file_pattern)
            form_pattern.envelope_from_shape(vault)
            form_pattern.selfweight_from_shape(vault)
            form_3D.selfweight_from_pattern(form_pattern)
        else:
            form_3D.selfweight_from_shape(vault)

        form_3D.envelope_from_shape(vault)

        # plot_form(form_3D, show_q=False).show()
        # form_3D.selfweight_from_shape(vault)

        # from compas_plotters import MeshPlotter
        # plotter = MeshPlotter(form_3D)
        # plotter.draw_edges()
        # plotter.draw_vertices(text={key: key for key in form_3D.vertices()})
        # plotter.show()

        if not arch_diagonal:
            if type_diagram_loads == 'fan_fd':
                keys_open_edge = [2, 3, 65, 125, 185, 245, 305, 365, 306, 246, 186, 126, 66, 5, 4]
            if type_diagram_loads == 'cross_fd':
                keys_open_edge = [28, 43, 58, 73, 88, 103, 118, 133, 148, 163, 178, 193, 208, 223, 224]
            if type_diagram_loads == 'topology-crossbraced':
                keys_open_edge = [1, 18, 35, 52, 69, 86, 103, 119, 134, 149, 164, 179, 194, 210, 224]
            loads[R] = [form_3D.vertex_attribute(key, 'pz') for key in keys_open_edge]
        if arch_diagonal:
            if type_diagram_loads == 'fan_fd':
                keys_diagonal = [2, 61, 121, 181, 241, 301, 361, 392, 363, 303, 243, 183, 123, 63, 6]
                crown = 392
            if type_diagram_loads == 'cross_fd':
                keys_diagonal = [28, 27, 41, 55, 69, 83, 97, 111, 125, 139, 153, 167, 181, 195, 209]
                crown = 111
            if type_diagram_loads == 'topology-crossbraced':
                break
            loads[R] = [form_3D.vertex_attribute(key, 'pz') if key != crown else 1/2*form_3D.vertex_attribute(key, 'pz') for key in keys_diagonal]

        # plotter = MeshPlotter(form_3D)
        # plotter.draw_edges()
        # plotter.draw_vertices(text={key: round(form_3D.vertex_attribute(key, 'pz'), 3) for key in form_3D.vertices()})
        # plotter.show()

        # from compas.datastructures import Mesh
        # for key in form_pattern.vertices():
        #     target = form_pattern.vertex_attribute(key, 'target')
        #     form_pattern.vertex_attribute(key, 'z', target)

        # dual = form_pattern.dual_diagram(Mesh)
        # # print(dual)
        # plotter = MeshPlotter(dual)
        # plotter.draw_edges()
        # plotter.draw_vertices(text={key: round(dual.vertex_attribute(key, 'z'), 2) for key in dual.vertices()})
        # plotter.show()

        # dual_fd = FormDiagram.from_mesh(dual)
        # for key in dual_fd.vertices():
        #     z = dual.vertex_attribute(key, 'z')
        #     dual_fd.vertex_attribute(key, 'z', z)
        # folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', 'pointed_crossvault', type_diagram_loads, 'tributary_area', 'R='+str(R), 'deg='+str(deg))
        # print(folder)
        # # os.mkdir(folder)
        # dual_address = os.path.join(folder, 'dual_' + 'pointed_crossvault' + '_' + type_diagram_loads + '_discr_' + str(discretisation) + '.json')
        # primal_address = os.path.join(folder, 'primal_' + 'pointed_crossvault' + '_' + type_diagram_loads + '_discr_' + str(discretisation) + '.json')
        # print(dual_address)
        # print(primal_address)
        # dual_fd.to_json(dual_address)
        # form_pattern.to_json(primal_address)
        # view_meshes([form_pattern, dual_fd]).show()

        # --------------------- 3.1 Create Minimisation for minimum thrust ---------------------

        optimiser = Optimiser()
        optimiser.data['library'] = 'Scipy'
        optimiser.data['solver'] = 'SLSQP'
        optimiser.data['constraints'] = ['funicular', 'envelope']  # 'reac_bounds'
        optimiser.data['variables'] = ['ind', 'zb', 't']
        optimiser.data['objective'] = 't'
        optimiser.data['printout'] = True
        optimiser.data['plot'] = False
        optimiser.data['find_inds'] = True
        optimiser.data['qmax'] = 5000.0
        optimiser.data['gradient'] = gradients
        optimiser.data['jacobian'] = gradients
        optimiser.data['thk'] = thk
        print(optimiser.data)

        # --------------------------- 3.2 Run optimisation with scipy ---------------------------

        analysis = Analysis.from_elements(arch, form, optimiser)
        # analysis.apply_selfweight()

        for index, load in enumerate(loads[R]):
            form.vertex_attribute(index, 'pz', load)

        # plotter = MeshPlotter(form)
        # plotter.draw_edges()
        # plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices()})
        # plotter.show()

        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()
        analysis.run()

        if optimiser.exitflag == 0:

            thk_min = form.attributes['thk']
            print(thk_min)
            sols[deg][R] = thk_min
            data_shape['thk'] = thk_min
            arch = Shape.from_library(data_shape)
            form.envelope_from_shape(arch)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', 'pointed_crossvault', type_diagram_loads, 'R='+str(R), 'min_thk')
            if deg:
                folder = os.path.join(folder, 'deg='+str(deg))
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            forms_address = os.path.join(folder, title)
            address_min = forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
            if save:
                form.to_json(address_min)

            # plot_reactions = 'black'  # or 'simple'
            # plot_form_xz(form, arch, show_q=False, fix_width=True, plot_reactions=plot_reactions, max_width=4,
            #             radius=0.09, stereotomy=blocks_on_plot, save=save_photo, cracks=False, hide_negative=True).show()

            save_photo = False
            blocks_on_plot = False
            if plot:
                plotter = plot_form_xz(form, arch, show_q=False, fix_width=True, plot_reactions=False, max_width=4,
                            radius=0.09, stereotomy=False, save=False, cracks=False, hide_negative=True)

                plotter.show()

            # view_solution(form).show()

print(sols)

for deg in sols:
    print(deg)
    for R in sols[deg]:
        print(R, ',', sols[deg][R])
