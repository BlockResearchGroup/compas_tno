import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
# from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
# from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints
from compas_tno.plotters import lookup_folder
from compas_tno.plotters import filter_min_thk
from compas_tno.plotters import plot_form_xz

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.75  # thickness on the start in meters
# Initial Settings
thk_reduction = 0.05  # in meters
thk_refined = 0.0001
limit_equal = 0.002
span = 10.0  # square span for analysis
k = 1
n = 1  # Discretisation for Surfaces...
Rs = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]
degs = [20]

# Rs = [6.5, 7.0]
# degs = [20]

# Rs = [8.5, 9.0, 9.5]
# degs = [20]

# FAN MINIMUMS
# degs = [0, 10, 20, 30, 40]
# Rs = [11.4027, 9.7364, 7.9422, 6.8896, 6.2453]

# files_fan = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=11.4027/min_thk/pointed_crossvault_fan_fd_discr_14_min_thk_32.59454210856643.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=9.7364/min_thk/deg=10/pointed_crossvault_fan_fd_discr_14_min_thk_28.040144250518782.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=7.9422/min_thk/deg=20/pointed_crossvault_fan_fd_discr_14_min_thk_21.06156962131167.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=6.8896/min_thk/deg=30/pointed_crossvault_fan_fd_discr_14_min_thk_15.38728690454782.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=6.2453/min_thk/deg=40/pointed_crossvault_fan_fd_discr_14_min_thk_11.060030746125904.json']

# CROSS MINIMUMS
degs = [0, 10, 20, 30, 40]
Rs = [7.2299, 6.7900, 6.1147, 5.6437, 5.3306]

# degs = [20]
# Rs = [6.1147]

# files_cross = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=7.2299/min_thk/pointed_crossvault_cross_fd_discr_14_min_thk_24.188374782556156.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.79/min_thk/deg=10/pointed_crossvault_cross_fd_discr_14_min_thk_20.339922896476374.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.1147/min_thk/deg=20/pointed_crossvault_cross_fd_discr_14_min_thk_13.878393453659262.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.6437/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_8.753882727013437.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.3306/min_thk/deg=40/pointed_crossvault_cross_fd_discr_14_min_thk_4.928791629011682.json']

minthk_list = []
addresses_main = []

# ---------------
type_formdiagram = 'cross_fd'
compare_with_arch = True
diagram_minimums = True
# ---------------

for j in range(len(degs)):
    addresses = []
    if diagram_minimums:
        radius_check = 1
    else:
        radius_check = len(Rs)

    for i in range(radius_check):
        if diagram_minimums:
            R = Rs[j]
        else:
            R = Rs[i]
        # R = Rs[i]
        # hc = hc_list[i]
        deg = degs[j]

        he = None

        sag = False
        smooth = False
        type_structure = 'pointed_crossvault'
        discretisation = 14

        radius = R
        A = span/(2*radius*(math.cos(math.radians(deg)) - 1) + span)
        xy_span_shape = [[-span/2*(A - 1), span*(1 + (A - 1)/2)], [-span/2*(A - 1), span*(1 + (A - 1)/2)]]
        hc_shape = A*math.sqrt(radius**2 - (radius - span/2)**2)

        if he:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'camb_h='+str(hc), 'min_thk')
        else:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(radius), 'min_thk')
        if deg:
            folder = os.path.join(folder, 'deg='+str(deg))
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
        if sag:
            title = title + 'sag_' + str(sag)
        if smooth:
            title = title + 'smooth_'
        forms_address = os.path.join(folder, title)

        print('**** Parameters [R, deg]:', R, deg)
        print('**** Folder:', folder)

        # data_shape = {
        #     'type': type_structure,
        #     'thk': thk,
        #     'discretisation': discretisation*n,
        #     'xy_span': xy_span_shape,
        #     't': 1.0,
        #     'hc': hc_shape,
        #     'hm': None,
        #     'he': he,
        # }

        files = lookup_folder(folder)
        # for f in files:
        #     os.rename(r'file path\OLD file name.file type',r'file path\NEW file name.file type')
        filtered = filter_min_thk(files, filters={'type_structure': 'pointed_crossvault', 'discretisation': discretisation, 'smooth': False, 'sag': False})
        limit_form_min = filtered[0]

        file_title = list(limit_form_min.keys())[0]
        thk_min = list(limit_form_min.values())[0]['thk']/100
        address = os.path.join(folder, file_title + '.json')

        print('**** Address:', address)
        print('**** Minimum thickness vault calculated:', thk_min, '\n')
        form = FormDiagram.from_json(address)
        minthk_list.append(thk_min)
        addresses.append(address)

        total_ind = 0
        for u, v in form.edges_where({'is_ind': True}):
            total_ind += 1
        print('Total Independent Edges:', total_ind)

        tol_cracks = 0.01 * thk_min

        plot_form(form, simple=True, show_q=False, cracks=True, tol_cracks=tol_cracks).show()

        # shape = Shape.from_formdiagram_and_attributes(form, data=data_shape)
        # view_solution(form).show()

        if compare_with_arch:

            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

            filtered = filter_min_thk(files, filters={'type_structure': 'pointed_arch', 'discretisation': discretisation})
            limit_form_min = filtered[0]

            file_title = list(limit_form_min.keys())[0]
            thk_min = list(limit_form_min.values())[0]['thk']/100
            address = os.path.join(folder, file_title + '.json')
            form_arch = FormDiagram.from_json(address)
            minthk_list.append(thk_min)
            addresses.append(address)

            print('**** Address:', address)
            print('**** Minimum thickness arch calculated:', thk_min, '\n')

            x0 = 0
            L = span
            t = 0.0
            b = 1.0
            A = L/(2*R*(math.cos(math.radians(deg)) - 1) + L)
            L_shape = A * L
            x0_shape = x0 + L/2*(1 - A)
            hc_shape = A*math.sqrt(R**2 - (R - L/2)**2)
            data_shape = {
                'type': 'pointed_arch',
                'hc': hc_shape,
                'L': L_shape,
                'thk': thk_min,
                'discretisation': discretisation + 1,
                'b': b,
                't': t,
                'x0': x0_shape
            }

            arch = Shape.from_library(data_shape)

            plot_form_xz(form, arch, show_q=False, fix_width=True, plot_reactions=False, max_width=4,
                radius=0.09, stereotomy=False, save=False, cracks=False, hide_negative=True, tol_cracks=tol_cracks).show()
            plot_form_xz(form_arch, arch, show_q=False, fix_width=True, plot_reactions=False, max_width=4,
                radius=0.09, stereotomy=False, save=False, cracks=False, hide_negative=True, tol_cracks=tol_cracks).show()

    addresses_main.append(addresses)

# print(minthk_list)
print(addresses_main)
