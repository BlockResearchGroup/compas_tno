import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
# from compas_tno.viewers import view_thrust
# from compas_tno.viewers import view_solution
# from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints
from compas_tno.plotters import lookup_folder
from compas_tno.plotters import filter_min_thk

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
Rs = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10]
hc_list = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
degs = [0, 10, 20, 30, 40]

addresses_main = []
minthk_list = []

# R = [5.0]
# hc_list = [5.0]

for j in range(len(degs)):
    addresses = []

    for i in range(len(Rs)):
        R = Rs[i]
        hc = hc_list[i]
        deg = degs[j]

        he = None

        sag = False
        smooth = False
        type_structure = 'pointed_crossvault'
        discretisation = 14
        type_formdiagram = 'fan_fd'

        radius = R
        A = span/(2*radius*(math.cos(math.radians(deg)) - 1) + span)
        xy_span_shape = [[-span/2*(A - 1), span*(1 + (A - 1)/2)], [-span/2*(A - 1), span*(1 + (A - 1)/2)]]
        hc_shape = A*math.sqrt(radius**2 - (radius - span/2)**2)

        if he:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'camb_h='+str(hc), 'min_thk')
        else:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc), 'min_thk')
        if deg:
            folder = os.path.join(folder, 'deg='+str(deg))
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
        if sag:
            title = title + 'sag_' + str(sag)
        if smooth:
            title = title + 'smooth_'
        forms_address = os.path.join(folder, title)

        print('**** Parameters [R, hc, deg]:', R, hc, deg)
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
        filtered = filter_min_thk(files, filters={'discretisation': discretisation, 'smooth': False, 'sag': False})
        limit_form_min = filtered[0]

        file_title = list(limit_form_min.keys())[0]
        thk_min = list(limit_form_min.values())[0]['thk']/100
        address = os.path.join(folder, file_title + '.json')

        print('**** Address:', address)
        print('**** Minimum thickness calculated:', thk_min, '\n')
        form = FormDiagram.from_json(address)
        minthk_list.append(thk_min)
        addresses.append(address)

        # shape = Shape.from_formdiagram_and_attributes(form, data=data_shape)
        # view_solution(form, shape).show()

    addresses_main.append(addresses)

# print(minthk_list)
print(addresses_main)
