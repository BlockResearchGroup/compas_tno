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

csv_file = os.path.join(folder, title + '_data.csv')
thicknesses, solutions = open_csv_row(csv_file)

# resave with different name
save_csv_row(thicknesses, solutions, path=csv_file, title=type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_sym')

img_graph = os.path.join(folder, title + '_diagram_symmetry.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True).show()

xy_limits = [[0.50, 0.15], [50, 0]]
img_graph = os.path.join(folder, title + '_diagram_limits_symmetry.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits).show()
