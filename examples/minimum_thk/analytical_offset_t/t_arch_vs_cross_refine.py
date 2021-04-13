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
import math

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN THICKNESS FOR POINTED ARCH --------------------
# ----------------------------------------------------------------------



sols = {}
loads = {}

Rs = [10.5, 12.0]
degs = [0, 10, 20, 30, 40]
discretisation = 14
thk = 0.75
b = 0.5  # Out of plane dimension  of arch
t = 0.0
type_diagram_loads = 'fan_fd'
type_structure = 'pointed_arch'
type_formdiagram = 'pointed_arch'

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

        data_diagram_3D = {
            'type': type_diagram_loads,
            'xy_span': [[0, L], [0, L]],
            'discretisation': discretisation,
            'fix': 'corners',
        }
        form_3D = FormDiagram.from_library(data_diagram_3D)
        # plot_form(form_3D, show_q=False).show()
        form_3D.selfweight_from_shape(vault)
        form_3D.envelope_from_shape(vault)
        keys_open_edge = [2, 3, 65, 125, 185, 245, 305, 365, 306, 246, 186, 126, 66, 5, 4]
        loads[R] = [form_3D.vertex_attribute(key, 'pz') for key in keys_open_edge]

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

        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()
        analysis.run()

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
        form.to_json(address_min)

        # plot_reactions = 'black'  # or 'simple'
        # plot_form_xz(form, arch, show_q=False, fix_width=True, plot_reactions=plot_reactions, max_width=4,
        #             radius=0.09, stereotomy=blocks_on_plot, save=save_photo, cracks=False, hide_negative=True).show()

        # plot = plot_form_xz(form, arch, show_q=False, fix_width=True, plot_reactions=False, max_width=4,
        #             radius=0.09, stereotomy=blocks_on_plot, save=save_photo, cracks=False, hide_negative=True)

        # plot.show()

print(sols)

for deg in sols:
    print(deg)
    for R in sols[deg]:
        print(R, ',', sols[deg][R])
