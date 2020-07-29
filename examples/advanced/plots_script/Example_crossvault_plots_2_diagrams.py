import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

xs = []
mins = []
maxs = []
legends = []

folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison')
title_bundled = 'Cross_vault_RQE'

type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Round Cross Vault - Cross_FD')

diagram_file = os.path.join(folder_main, title_bundled + '_diagramA.pdf')
colors = ['#D2B4DE', '#b4c0de'][:len(legends)]
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.15, 0.04], [130, 50]]).show()


type_structure = 'crossvault'
type_formdiagram = 'fan_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Round Cross Vault - Fan_FD')

diagram_file = os.path.join(folder_main, title_bundled + '_diagramB.pdf')
colors = ['#D2B4DE', '#b4c0de'][:len(legends)]
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.15, 0.04], [130, 50]]).show()
