from compas_tno.plotters import open_csv
import os
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import save_csv_row
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_thrust

type_structure = 'pavillionvault'
type_formdiagram = 'ortho_fd'  # Try also 'fan_fd'
discretisation = 10
hc = 5.0

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# ----------------------- Save output data --------------------------

# csv_file = os.path.join(folder, title + '_data.csv')

# csv_file = '/Users/mricardo/compas_dev/me/minmax/2D_arch/diagram_thrust/diagram_thrust_arch.csv'
# thicknesses, solutions = open_csv_row(csv_file)
# x, min_thu, max_thu = open_csv(csv_file)
# print(thicknesses)

# print(solutions)

thk = [0.2, 0.195, 0.19, 0.185, 0.18, 0.175, 0.17, 0.165, 0.16, 0.155, 0.15, 0.145, 0.14, 0.135, 0.13, 0.125, 0.12, 0.115, 0.11, 0.1079]
min_thu = [0.316, 0.319, 0.323, 0.327, 0.331, 0.335, 0.339, 0.343, 0.347, 0.351, 0.356, 0.36, 0.365, 0.369, 0.374, 0.379, 0.384, 0.389, 0.394, 0.396]
max_thu = [0.512, 0.507, 0.502, 0.496, 0.491, 0.485, 0.48, 0.474, 0.469, 0.462, 0.455, 0.449, 0.442, 0.436, 0.43, 0.423, 0.417, 0.408, 0.4, 0.397]

# for i in range(len(max_thu)):
#     max_thu[i] *= -1

thicknesses = [thk, thk]
solutions = [max_thu, min_thu]

# resave with different name
# save_csv_row(thicknesses, solutions, path=csv_file, title=type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_sym')

img_graph = None # os.path.join(folder, title + '_diagram_symmetry.pdf')
# diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True).show()

xy_limits = [[0.20, 0.10], [60, 30]]
# img_graph = os.path.join(folder, title + '_diagram_limits_symmetry.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, GSF_ticks=[1.0, 1.5, 1.75, 2.0], xy_limits=xy_limits).show()
