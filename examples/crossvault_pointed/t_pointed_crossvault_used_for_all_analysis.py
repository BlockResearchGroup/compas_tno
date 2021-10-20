import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
from compas_tno.viewers import view_thrust
from compas_tno.viewers import Viewer
from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints
from compas_tno.plotters import lookup_folder
from compas_tno.plotters import filter_min_thk
import time
os.environ['QT_MAC_WANTS_LAYER'] = '1'

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THK FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
# Initial Settings
thk_reduction = 0.05  # in meters
thk_refined = 0.0001
limit_equal = 0.002
minimise_thickness = True
discretisation = 20

span = 10.0  # square span for analysis
k = 1
n = 1  # Discretisation for Surfaces...

# degs = [40]
# R = [5.3306]

# R = [5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]
# degs = [20]

# R = [7.080]
# degs = [20]

# cross minimum
# degs = [20]
# R = [6.1147]

# fan minimum
# degs = [20]
# R = [7.9422]

# Paper study
degs = [20]
R = [6.1147, 7.08, 7.9422]

# ---------------------------------------
t0 = thk = 0.50  # thickness on the start in meters
he = None  # [5.0, 5.0, 5.0, 5.0]
sols = {}
time_calc = {}
sag = False  # 10.0 50.0
smooth = False
lookup_refine = False
save = True
compute_diagram_of_thrust = False
first_plot = False
type_formdiagram_pattern = False #'topology-crosssmooth2'
solver = 'SLSQP'
diagrams = ['cross_fd']
# ---------------------------------------

# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

for type_formdiagram in diagrams:

    sols[type_formdiagram] = {}
    time_calc[type_formdiagram] = {}

    for deg in degs:

        type_structure = 'pointed_crossvault'
        # type_formdiagram = 'cross_fd'  # Try also 'fan_fd
        gradients = True

        sols[type_formdiagram][deg] = {}
        time_calc[type_formdiagram][deg] = {}

        if '-' in [char for char in type_formdiagram]:
            type_topology = type_formdiagram.split('-')[-1]
        else:
            type_topology = False

        if type_formdiagram == 'topology-mix':
            smooth = True

        if type_topology:
            type_formdiagram = 'topology-' + type_topology
            fd_mesh = 'FormDiagram-' + type_topology

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
            if first_plot:
                plot_form(form, show_q=False, fix_width=10).show()
            apply_sag(form, boundary_force=sag)
        if smooth:
            print('Smoothing ON')
            if first_plot:
                plot_form(form, show_q=False, fix_width=10).show()
            cons = {}
            for key in form.vertices_on_boundary():
                cons[key] = form.vertex_coordinates(key)
            # for key in form.vertices():
            #     x, y, z = form.vertex_coordinates(key)
            #     if abs(x - 5.0) < 0.01 or abs(y - 5.0) < 0.01:
            #         cons[key] = form.vertex_coordinates(key)
            # cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
            constrained_smoothing(form, damping=0.5, kmax=300, constraints=cons, algorithm='centroid')

        if first_plot:
            plot_form(form, show_q=False, fix_width=10).show()
            first_plot = False

        # --------------------- Create Optimiser ---------------------

        optimiser = Optimiser()
        if solver == 'SLSQP':
            optimiser.settings['library'] = 'Scipy'
            optimiser.settings['solver'] = 'SLSQP'
        else:
            optimiser.settings['library'] = 'IPOPT'
            optimiser.settings['solver'] = 'IPOPT'
        optimiser.settings['max_iter'] = 2000
        optimiser.settings['constraints'] = ['funicular', 'envelope']
        optimiser.settings['printout'] = True  # True
        optimiser.settings['plot'] = False
        optimiser.settings['find_inds'] = True
        optimiser.settings['qmax'] = 10e+10
        optimiser.settings['gradient'] = gradients
        optimiser.settings['jacobian'] = gradients
        optimiser.settings['max_thk'] = thk * 1.0
        if minimise_thickness:
            optimiser.settings['objective'] = 't'
            optimiser.settings['variables'] = ['ind', 'zb', 't']
        else:
            optimiser.settings['objective'] = 'bestfit'
            optimiser.settings['variables'] = ['ind', 'zb']

        for i in range(len(R)):

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
                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'camb_h='+str(hc_shape), 'min_thk')
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

            if lookup_refine:

                type_formdiagram_load = type_formdiagram
                folder_pattern = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram_load, 'R='+str(radius), 'min_thk')
                if deg:
                    folder_pattern = os.path.join(folder_pattern, 'deg='+str(deg))

                print('Lookup Folder:', folder_pattern)
                files = lookup_folder(folder_pattern)
                filtered = filter_min_thk(files, filters={'type_structure': 'pointed_crossvault', 'objective': 'min', 'discretisation': discretisation})
                limit_form_min = filtered[0]

                file_title = list(limit_form_min.keys())[0]
                thk_min = list(limit_form_min.values())[0]['thk']/100
                address = os.path.join(folder, file_title + '.json')
                # address_pattern = os.path.join(folder_pattern, file_title + '.json')

                print('**** Address:', address)
                # print('**** Address-Pattern:', address_pattern)
                print('**** Minimum thickness vault Filtered:', thk_min, '\n')
                form = FormDiagram.from_json(address)
                # form_pattern = FormDiagram.from_json(address_pattern)
                # form_pattern.selfweight_from_shape(vault)
                # form.selfweight_from_pattern(form_pattern)
                # form.envelope_from_shape(vault)

                compute_initial_point = False

                optimiser.settings['max_thk'] = thk_min
                # if radius == 7.5:
                #     optimiser.settings['max_thk'] = 0.2419
                #     optimiser.settings['max_iter'] = 500
                tmin_reduced = thk_min
                print('Constraining it for:', tmin_reduced)
                thk = tmin_reduced

                # plot_form(form, show_q=False, cracks=True).show()
                # view_thrust(form).show()
                # view_solution(form, vault).show()

            else:

                compute_initial_point = True

            # --------------------- Create Initial point with TNA ---------------------

            if compute_initial_point:

                address_lp = forms_address + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'
                address_shape = forms_address + '_' + 'shape' + '_thk_' + str(100*thk) + '.json'

                form.selfweight_from_shape(vault)
                form.envelope_from_shape(vault)
                form.to_json(address_shape)
                print('\n', address_shape, '\n')

                break
                s

                if type_formdiagram_pattern:
                    file_pattern = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram_pattern, 'FormDiagram-' + type_formdiagram_pattern.split('-')[-1] + '.json')
                    form_pattern = FormDiagram.from_json(file_pattern)
                    form_pattern.selfweight_from_shape(vault)
                    form.selfweight_from_pattern(form_pattern)
                else:
                    form.selfweight_from_shape(vault)

                # from compas_plotters import MeshPlotter
                # plotter = MeshPlotter(form_pattern_loads, figsize=(10, 10))
                # plotter.draw_edges()
                # plotter.draw_vertices(text={key: round(form_pattern_loads.vertex_attribute(key, 'pz'), 3) for key in form_pattern_loads.vertices()})
                # # plotter.draw_vertices(text={key: key for key in form.vertices()})
                # # plotter.draw_faces(text={key: key for key in form.faces()})
                # plotter.show()

                # from compas_plotters import MeshPlotter
                # plotter = MeshPlotter(form, figsize=(10, 10))
                # plotter.draw_edges()
                # plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 3) for key in form.vertices()})
                # # plotter.draw_vertices(text={key: key for key in form.vertices()})
                # # plotter.draw_faces(text={key: key for key in form.faces()})
                # plotter.show()

                try:
                    form = FormDiagram.from_json(address_lp)
                    print('Loaded LP form diagram')
                except BaseException:
                    print('LP form diagram not found')
                    print(address_lp)
                    # form.selfweight_from_shape(vault)
                    form.envelope_from_shape(vault)
                    form.initialise_loadpath()
                    form.to_json(address_lp)

                # starting_point_file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=6.0/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_17.723600813371625.json'
                # starting_point_file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=8.0/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_27.437516049296605.json'
                # starting_point_file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=7.9422/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_27.047244399550735.json'
                # form = FormDiagram.from_json(starting_point_file)

                # file_pattern = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', 'pointed_crossvault', 'topology-crosssmooth2', 'FormDiagram-crosssmooth2.json')
                # form_pattern = FormDiagram.from_json(file_pattern)

                # form.form_update_with_parallelisation()

                # plot_form(form, show_q=False).show()

                # from compas_tno.algorithms import z_from_form
                # z_from_form(form)

                # plot_form(form, show_q=False).show()
                # view_thrust(form).show()
                # view_solution(form, vault).show()

            if compute_diagram_of_thrust:

                # ---------------------- Analysis for the Diagram of Thrust -------------

                analysis = Analysis.from_elements(vault, form, optimiser)
                # results = analysis.thk_minmax_GSF(thk, thk_step=thk_reduction, save_forms=forms_address,
                #                 plot=False, printout=False, jump_minthk=True, swt_from_pattern=form_pattern)
                results = analysis.thk_minmax_GSF(thk, thk_step=thk_reduction, save_forms=forms_address,
                                plot=False, printout=False)
                thicknesses, solutions = results

                # ----------------------- Save output data --------------------------

                csv_file = os.path.join(folder, title + '_data.csv')
                save_csv_row(thicknesses, solutions, path=csv_file, title=title)
                print('Saved csv file:', csv_file)

                img_graph = os.path.join(folder, title + '_diagram.pdf')
                diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True)

                xy_limits = [[0.60, 0.20], [120, 30]]
                img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
                diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits)

                thk_min = thicknesses[0][-1]
                sols[type_formdiagram][deg][radius] = thk_min

            else:

                analysis = Analysis.from_elements(vault, form, optimiser)
                # analysis.apply_selfweight()
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
                    address_min_check = forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
                    # address_max = forms_address + '_' + 'max' + '_thk_' + str(100*thk_min) + '.json'

                    if optimiser.exitflag == 0:
                        # form.to_json(address_max)
                        if save:
                            form.to_json(address_min)
                            print('Salved form to:', address_min)
                        else:
                            print('SOLVED, but did not save to:', address_min)
                        sols[type_formdiagram][deg][radius] = thk_min
                        time_calc[type_formdiagram][deg][radius] = optimiser.time
                    else:
                        form.to_json(address_min_check)
                        print('Salved TO BE CHECKED FORM TO:', address_min_check)
                        sols[type_formdiagram][deg][radius] = [thk_min]
                        time_calc[type_formdiagram][deg][radius] = optimiser.time

                    plot_form(form, simple=True, show_q=False, cracks=True).show()
                    # plot_form(form, simple=True, show_q=False, cracks=True, save=forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.pdf').show()
                    # view_thrust(form).show()
                    # view = Viewer(form)
view.show_solution()

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


# print('******** Resume of Times')
# for key in time_calc:
#     print(key)
#     for discr in time_calc[key]:
#         print(discr)
#         for hc in time_calc[key][discr]:
#             print(hc, ',', time_calc[key][discr][hc])

