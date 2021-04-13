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
from compas_tno.plotters import lookup_folder
from compas_tno.plotters import filter_min_thk
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.20  # thickness on the start in meters
# Initial Settings
thk_reduction = 0.05  # in meters
thk_refined = 0.0001
limit_equal = 0.002
minimise_thickness = True
discretisation = 14

refine = True
refine_limit = 0.0005

span = 10.0  # square span for analysis
k = 1
n = 1  # Discretisation for Surfaces...
R = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10]
hc_list = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
degs = [0, 10, 20, 30, 40]

# degs = [0]
# R = [6.5, 8.0]

# degs = [10]
# R = [6.0, 7.0]

# degs = [20]
# R = [5.5, 6.5]

# degs = [30]
# R = [5.0, 6.0]

degs = [20]
R = [6.0, 7.0]

# 10 - 7 to 8
# 20 - 6 to 7
# 30 - 5.5 to 6.5
# 40 - 5 to 6
min_thk = [0, 0]
# hc_list = [5.0]

golden = (math.sqrt(5) - 1)/2

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
type_formdiagram_pattern = 'topology-crosssmooth2'

# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

for type_formdiagram in ['topology-crossbraced']:  # THIS CAN ONLY WORK IF WE ADD THAT A DIFFERENT SWT MUST BE CONSIDERED ;)

    sols[type_formdiagram] = {}

    for deg in degs:

        type_structure = 'pointed_crossvault'
        # type_formdiagram = 'cross_fd'  # Try also 'fan_fd
        gradients = True

        sols[type_formdiagram][deg] = {}

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
            apply_sag(form, boundary_force=sag)
        if smooth:
            print('Smoothing ON')
            cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
            constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

        # plot_form(form, show_q=False, fix_width=10).show()

        # --------------------- Create Optimiser ---------------------

        optimiser = Optimiser()
        # optimiser.data['library'] = 'Scipy'
        # optimiser.data['solver'] = 'SLSQP'
        optimiser.data['library'] = 'IPOPT'
        optimiser.data['solver'] = 'IPOPT'
        optimiser.data['max_iter'] = 1000
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

        j = 0
        for radius in R:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(radius), 'min_thk')
            if deg:
                folder = os.path.join(folder, 'deg='+str(deg))
            print('------ Looking up folder:', folder)
            files = lookup_folder(folder)
            if type_formdiagram == 'topology-mix':
                filtered = filter_min_thk(files, filters={'type_structure': 'pointed_crossvault', 'smooth': smooth, 'sag': False})
            else:
                filtered = filter_min_thk(files, filters={'type_structure': 'pointed_crossvault', 'discretisation': discretisation, 'smooth': smooth, 'sag': False})
            limit_form_min = filtered[0]
            file_title = list(limit_form_min.keys())[0]
            thk_min = list(limit_form_min.values())[0]['thk']/100
            address = os.path.join(folder, file_title + '.json')
            print('**** Address:', address)
            print('**** Minimum thickness calculated:', thk_min, '\n')
            form = FormDiagram.from_json(address)
            min_thk[j] = thk_min
            sols[type_formdiagram][deg][radius] = thk_min
            j += 1

        print(sols)

        fx1 = None
        fx2 = None
        x1 = R[0] + (R[1] - R[0]) * golden
        x2 = R[1] - (R[1] - R[0]) * golden
        i = 0

        print('Begining process of this iteration / R / x1 / x2:', R, x1, x2)
        print('Values on Beginning / min thk / fx1 / fx2:', min_thk, fx1, fx2)

        while refine:

            print('\n\n---------------')
            print('Iteration i:', i)

            if i == 0:
                radius_to_check = [x1, x2]
            else:
                if fx1 > fx2:   # vanish b | x1 become x2 | x2 calculated
                    print('fx1 > fx2', fx1, fx2)
                    print('x1, x2', x1, x2)
                    R[1] = x1
                    min_thk[1] = fx1
                    x1 = x2
                    fx1 = fx2
                    x2 = R[1] - (R[1] - R[0]) * golden
                    fx2 = None
                    radius_to_check = [x2]
                else:           # vanish a | x2 become x1 | x1 calculated
                    print('fx1 < fx2', fx1, fx2)
                    print('x1, x2', x1, x2)
                    R[0] = x2
                    min_thk[0] = fx2
                    x2 = x1
                    fx2 = fx1
                    x1 = R[0] + (R[1] - R[0]) * golden
                    fx1 = None
                    radius_to_check = [x1]

            print('Bounds of the Problem:', R)
            print('Radius to Check:', radius_to_check)

            if R[1] - R[0] < refine_limit:
                print('Breaking due to limit', R, refine_limit)
                refine = False
                radius_to_check = [(R[0] + R[1] + x1 + x2)/4]

            # if i > 20:
            #     print('Breaking due maximum i:', R, i)
            #     refine = False
            #     radius_to_check = [(R[0] + R[1] + x1 + x2)/4]

            for radius in radius_to_check:

                A = span/(2*radius*(math.cos(math.radians(deg)) - 1) + span)
                xy_span_shape = [[-span/2*(A - 1), span*(1 + (A - 1)/2)], [-span/2*(A - 1), span*(1 + (A - 1)/2)]]
                hc_shape = A*math.sqrt(radius**2 - (radius - span/2)**2)
                if he:
                    he = [A*he[0], A*he[1], A*he[2], A*he[3]]

                print('\n**** Starting Analysis for: radius={0} - thk={1} - deg={2} ****\n'.format(radius, thk, deg))

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

                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'refined', 'min_thk')
                if deg:
                    folder = os.path.join(folder, 'deg='+str(deg))
                os.makedirs(folder, exist_ok=True)
                title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_R=' + str(radius)
                if sag:
                    title = title + 'sag_' + str(sag)
                if smooth:
                    title = title + 'smooth_'
                forms_address = os.path.join(folder, title)

                # --------------------- Create Initial point with TNA ---------------------

                if i == 0:
                    address_lp = forms_address + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'

                    if type_formdiagram_pattern:
                        file_pattern = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram_pattern, 'FormDiagram-' + type_formdiagram_pattern.split('-')[-1] + '.json')
                        form_pattern = FormDiagram.from_json(file_pattern)
                        form_pattern.selfweight_from_shape(vault)
                        form.selfweight_from_pattern(form_pattern, plot=True)
                    else:
                        form.selfweight_from_shape(vault)
                    form.envelope_from_shape(vault)

                    try:
                        form = FormDiagram.from_json(address_lp)
                        print('Loaded LP form diagram')
                    except:
                        print('LP form diagram not found')
                        print(address_lp)
                        # form.selfweight_from_shape(vault)
                        # form.envelope_from_shape(vault)
                        # form = form.initialise_tna(plot=False)
                        form.initialise_loadpath()
                        # plot_form(form, show_q=False).show()
                        address_lp = forms_address + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'
                        form.to_json(address_lp)

                analysis = Analysis.from_elements(vault, form, optimiser)
                # analysis.apply_selfweight()
                analysis.apply_envelope()
                analysis.apply_reaction_bounds()
                analysis.set_up_optimiser()
                analysis.run()

                thk_min = form.attributes['thk']
                data_shape['thk'] = thk_min
                print('Minimum thickness calculated = ', thk_min)
                vault = Shape.from_library(data_shape)
                form.envelope_from_shape(vault)

                address_min = forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
                # address_max = forms_address + '_' + 'max' + '_thk_' + str(100*thk_min) + '.json'

                if optimiser.exitflag == 0:
                    # form.to_json(address_max)
                    if not refine:
                        form.to_json(address_min)

                sols[type_formdiagram][deg][radius] = thk_min

                if not refine:
                    plot_form(form, simple=True, show_q=False, cracks=True, save=forms_address + '_' + 'min' + '_thk_' + str(100*thk_min) + '.png').show()
                # view_thrust(form).show()
                # view_solution(form, vault).show()

                if i == 0:
                    if not fx1:
                        fx1 = thk_min
                    elif not fx2:
                        fx2 = thk_min
                    else:
                        raise Exception
                else:
                    if not fx1:
                        fx1 = thk_min
                    if not fx2:
                        fx2 = thk_min

            print('\nEnd of this iteration / R / x1 / x2:', R, x1, x2)
            print('Values of this iteration / min thk / fx1 / fx2:', min_thk, fx1, fx2)

            i = i + 1

            print(sols)

print('******** Resume of Solutions')
for key in sols:
    print(key)
    for discr in sols[key]:
        print(discr)
        for hc in sols[key][discr]:
            print(hc, ',', sols[key][discr][hc])
