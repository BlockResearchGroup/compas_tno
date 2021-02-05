import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
# Initial Settings
thk_reduction = 0.05  # in meters
thk_refined = 0.0001
limit_equal = 0.002
minimise_thickness = True
discretisation = 14

span = 10.0  # square span for analysis
k = 1
n = 1  # Discretisation for Surfaces...
R = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10]
hc_list = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
degs = [0, 10, 20, 30, 40]

degs = [30]
R = [6.8692]
# hc_list = [5.0]

# i =8
# R = [R[i]]
# hc_list = [hc_list[i]]
# degs = [40]

# hc_list = [6.71]
he = None  # [5.0, 5.0, 5.0, 5.0]
sols = {}
# [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
sag = False
smooth = False

compute_diagram_of_thrust = False

# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

for type_formdiagram in ['fan_fd']:

    sols[type_formdiagram] = {}

    for deg in degs:

        type_structure = 'pointed_crossvault'
        # type_formdiagram = 'cross_fd'  # Try also 'fan_fd
        gradients = True

        sols[type_formdiagram][deg] = {}

        if type_formdiagram == 'topology-mix':
            type_topology = 'mix'  # 'mix' or None
        else:
            type_topology = None
        if type_topology:
            type_formdiagram = 'topology-' + type_topology
            fd_mesh = 'FormDiagram-' + type_topology
            smooth = True

        # ----------------------- Create Form Diagram for analysis ---------------------------

        if not type_topology:
            data_diagram = {
                'type': type_formdiagram,
                'xy_span': [[0, span], [0, k*span]],
                'discretisation': discretisation,
                'fix': 'corners',
            }
            form = FormDiagram.from_library(data_diagram)
        else:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
            form = FormDiagram.from_json((os.path.join(folder, fd_mesh + '.json')))

        if sag:
            print('Sag ON')
            apply_sag(form, boundary_force=sag)
        if smooth:
            print('Smoothing ON')
            cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
            constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

        # plot_form(form, show_q=False, fix_width=10).show()

        # --------------------- Create Optimiser ---------------------

        optimiser = Optimiser()
        optimiser.data['library'] = 'Scipy'
        optimiser.data['solver'] = 'SLSQP'
        # optimiser.data['library'] = 'IPOPT'
        # optimiser.data['solver'] = 'IPOPT'
        optimiser.data['constraints'] = ['funicular', 'envelope']
        optimiser.data['printout'] = True
        optimiser.data['plot'] = False
        optimiser.data['find_inds'] = True
        optimiser.data['qmax'] = 10e+10
        optimiser.data['gradient'] = gradients
        optimiser.data['jacobian'] = gradients
        if minimise_thickness:
            optimiser.data['objective'] = 't'
            optimiser.data['variables'] = ['ind', 'zb', 't']
        else:
            optimiser.data['objective'] = 'bestfit'
            optimiser.data['variables'] = ['ind', 'zb']

        for i in range(len(R)):

            # hc = hc_list[i]
            radius = R[i]
            A = span/(2*radius*(math.cos(math.radians(deg)) - 1) + span)
            xy_span_shape = [[-span/2*(A - 1), span*(1 + (A - 1)/2)], [-span/2*(A - 1), span*(1 + (A - 1)/2)]]
            hc_shape = A*math.sqrt(radius**2 - (radius - span/2)**2)
            if he:
                he = [A*he[0], A*he[1], A*he[2], A*he[3]]

            print('\n**** Starting Analysis for: radius={0} and thk={1} ****\n'.format(radius, thk))
            if deg:
                print('---- Considering deg ={0} ****\n'.format(deg))

            # --------------------- Shape with initial THK ---------------------

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'discretisation': discretisation,
                'xy_span': xy_span_shape,
                't': 0.0,
                'hc': hc_shape,
                'hm': None,
                'he': he,
            }

            vault = Shape.from_library(data_shape)
            # view_shapes(vault).show()

            # ----------------------- Create Analysis loop on limit analysis --------------------------

            if he:
                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'camb_h='+str(hc), 'min_thk')
            else:
                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(radius), 'min_thk')
            if deg:
                folder = os.path.join(folder, 'deg='+str(deg))
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            if sag:
                title = title + 'sag_' + str(sag)
            if smooth:
                title = title + 'smooth_'
            forms_address = os.path.join(folder, title)

            # --------------------- Create Initial point with TNA ---------------------

            address_lp = forms_address + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'
            # address_shape = forms_address + '_' + 'shape' + '_thk_' + str(100*thk) + '.json'

            form.selfweight_from_shape(vault)
            form.envelope_from_shape(vault)
            # form.to_json(address_shape)
            # print('\n', address_shape, '\n')

            try:
                form = FormDiagram.from_json(address_lp)
                print('Loaded LP form diagram')
            except:
                print('LP form diagram not found')
                print(address_lp)
                form.selfweight_from_shape(vault)
                form.envelope_from_shape(vault)
                # form = form.initialise_tna(plot=False)
                form.initialise_loadpath()
                # plot_form(form, show_q=False).show()
                address_lp = forms_address + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'
                form.to_json(address_lp)

            # view_thrust(form).show()
            # view_solution(form, vault).show()

            if compute_diagram_of_thrust:

                # ---------------------- Analysis for the Diagram of Thrust -------------

                analysis = Analysis.from_elements(vault, form, optimiser)
                results = analysis.thk_minmax_GSF(thk, thk_step=thk_reduction, save_forms=forms_address, plot=False, printout=False)
                thicknesses, solutions = results

                # ----------------------- Save output data --------------------------

                csv_file = os.path.join(folder, title + '_data.csv')
                save_csv_row(thicknesses, solutions, path=csv_file, title=title)

                img_graph = os.path.join(folder, title + '_diagram.pdf')
                diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True)

                xy_limits = [[0.60, 0.20], [120, 30]]
                img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
                diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits)

                thk_min = thicknesses[0][-1]
                sols[type_formdiagram][deg][hc] = thk_min

            else:

                analysis = Analysis.from_elements(vault, form, optimiser)
                analysis.apply_selfweight()
                analysis.apply_envelope()
                analysis.apply_reaction_bounds()
                analysis.set_up_optimiser()
                analysis.run()

                if minimise_thickness:

                    thk_min = form.attributes['thk']
                    data_shape['thk'] = thk_min
                    print('Minimum thickness calculated = ', thk_min)
                    vault = Shape.from_library(data_shape)
                    form.envelope_from_shape(vault)

                    address_min = forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
                    # address_max = forms_address + '_' + 'max' + '_thk_' + str(100*thk_min) + '.json'

                    if optimiser.exitflag == 0:
                        # form.to_json(address_max)
                        form.to_json(address_min)

                    sols[type_formdiagram][deg][radius] = thk_min

                    plot_form(form, simple=True, show_q=False, cracks=True, save=forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.pdf').show()
                    # view_thrust(form).show()
                    # view_solution(form, vault).show()

                else:

                    distance_target = optimiser.fopt  # form.attributes['fopt']
                    print('Distance Squared = ', distance_target)

                    address_bestfit = forms_address + '_' + 'bestfit' + '_thk_' + str(100*thk) + '.json'

                    if optimiser.exitflag == 0:
                        form.to_json(address_bestfit)

                    sols[type_formdiagram][deg][radius] = distance_target

        print(sols)

print('******** Resume of Solutions')
for key in sols:
    print(key)
    for discr in sols[key]:
        print(discr)
        for hc in sols[key][discr]:
            print(hc, ',', sols[key][discr][hc])