import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

xs = []
mins = []
maxs = []
legends = []

type_structure = 'crossvault'
type_formdiagram = 'cross_diagonal'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

type_structure = 'crossvault'
type_formdiagram = 'fan_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_diagonal'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'
discretisation = 10
folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append(type_structure + '\n' + type_formdiagram)

# type_structure = 'dome'
# type_formdiagram = 'radial_fd'
# discretisation = [10,20]
# folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# csv_file = os.path.join(folder, title + '_data.csv')
# print(csv_file)
# size_parameters, solutions_min, solutions_max = open_csv(csv_file)

# size_parameters = [0.15, 0.13999999999999999, 0.13, 0.12, 0.11000000000000001, 0.1, 0.09, 0.08, 0.06999999999999999, 0.06, 0.05, 0.04]
# solutions_min = [0.17731215042243287, 0.17806468465527767, 0.18305438897457416, 0.18904946660086341, 0.19516219439387025, 0.2014276343410696, 0.20787747820215766, 0.2145417988278732, 0.22145024547358474, 0.22927447092577885, 0.23762113394134857, 0.24641356394496203]
# solutions_max = [-0.7155566816863638, -0.630978436849197, -0.5636444115639517, -0.5080959420004013, -0.46103526483892404, -0.42034510711099815, -0.38458914105216374, -0.35275371312098464, -0.3240995465603163, -0.2980725812098216, -0.27424792774876355, -0.2522933017928862]

# xs.append(size_parameters)
# mins.append(solutions_min)
# maxs.append(solutions_max)
# legends.append('Hemispherical Dome')

folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data')
title_bundled = 'NewFDs'
diagram_file = os.path.join(folder, title_bundled + '_diagram.pdf')
print(diagram_file)
colors = ['#17202A', '#D2B4DE', '#2ECC71', '#17202A', '#D2B4DE', '#2ECC71']
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, xy_limits=[[0.10, 0.04], [110, 40]]).show()
