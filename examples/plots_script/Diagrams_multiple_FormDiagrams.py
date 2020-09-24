import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv
from compas_tno.plotters import open_csv_row

xs = []
mins = []
maxs = []
legends = []

folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison')
title_bundled = 'All_Structures'

type_structures = ['pointed_crossvault']*4
type_formdiagrams = ['cross_fd', 'topology-mix', 'fan_fd',]
discretisations = [10]*4
# colors = ['#D2B4DE', '#F39C12', '#2ECC71', '#17202A', '#c7b66f']
colors = ['#F39C12', '#D2B4DE', '#2ECC71', '#17202A', '#c7b66f']
legends = []
xy_limits = [[0.60, 0.05], [130, 10]]

thicknesses_combined = []
solutions_combined = []

for i in range(len(type_formdiagrams)):
    type_structure = type_structures[i]
    type_formdiagram = type_formdiagrams[i]
    discretisation = discretisations[i]
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h=5.0')
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    if type_formdiagram == 'topology-mix':
        title = 'pointed_crossvault_topology-mix_Mesh-mix_discr_100smooth_'
    csv_file = os.path.join(folder, title + '_data.csv')
    print(csv_file)
    thicknesses, solutions = open_csv_row(csv_file)
    thicknesses_combined.append(thicknesses)
    solutions_combined.append(solutions)
    legends.append(type_structure + '-' + type_formdiagram)
    diagram_file = os.path.join(folder_main, str(legends) + '.pdf')
    diagram_of_multiple_thrust(thicknesses_combined, solutions_combined, legends, colors=colors, save=diagram_file, fill=True, show_legend=True, xy_limits=xy_limits).show()


# type_structure = 'pointed_crossvault'
# type_formdiagram = 'cross_fd'
# discretisation = 10
# folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# csv_file = os.path.join(folder, title + '_data.csv')
# print(csv_file)
# size_parameters, solutions_min, solutions_max = open_csv(csv_file)
# xs.append(size_parameters)
# mins.append(solutions_min)
# maxs.append(solutions_max)
# legends.append('Gothic Cross Vault')

# diagram_file = os.path.join(folder_main, title_bundled + '_diagramC.pdf')
# colors = ['#D2B4DE', '#F39C12', '#2ECC71', '#17202A'][:len(legends)]
# diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.15, 0.02], [130, 10]]).show()

# # type_structure = 'dome'
# # type_formdiagram = 'radial_fd'
# # discretisation = [10,20]
# # folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
# # title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# # csv_file = os.path.join(folder, title + '_data.csv')
# # print(csv_file)
# # size_parameters, solutions_min, solutions_max = open_csv(csv_file)

# # dome_radial_spaced_fd_discr_[8, 20]_data.csv

# size_parameters = [0.15, 0.13999999999999999, 0.13, 0.12, 0.11000000000000001, 0.1, 0.09, 0.08, 0.06999999999999999, 0.06, 0.05, 0.04]
# solutions_min = [0.17731215042243287, 0.17806468465527767, 0.18305438897457416, 0.18904946660086341, 0.19516219439387025, 0.2014276343410696, 0.20787747820215766, 0.2145417988278732, 0.22145024547358474, 0.22927447092577885, 0.23762113394134857, 0.24641356394496203]
# solutions_max = [-0.7155566816863638, -0.630978436849197, -0.5636444115639517, -0.5080959420004013, -0.46103526483892404, -0.42034510711099815, -0.38458914105216374, -0.35275371312098464, -0.3240995465603163, -0.2980725812098216, -0.27424792774876355, -0.2522933017928862]


# # size_parameters = [0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05]
# # solutions_min = [0.171, 0.172, 0.176, 0.182, 0.189, 0.195, 0.202, 0.210, 0.217, 0.224, 0.232]
# # solutions_max = [-2.478, -2.655, -2.860, -3.098, -3.186, -3.717, -1.647, -0.876, -0.563, -0.415, -0.314]

# # 0.171, 0.172, 0.176, 0.182, 0.189, 0.195, 0.202, 0.210, 0.217, 0.224, 0.232
# # -2.478, -2.655, -2.860, -3.098,-3.186, -3.717, -1.647, -0.876, -0.563, -0.415, -0.314

# xs.append(size_parameters)
# mins.append(solutions_min)
# maxs.append(solutions_max)

# legends.append('Hemispherical Dome')

# diagram_file = os.path.join(folder_main, title_bundled + '_diagramD.pdf')
# print(diagram_file)
# colors = ['#D2B4DE', '#F39C12', '#2ECC71', '#17202A']
# diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.15, 0.02], [130, 10]]).show()


# type_structure = 'arch'
# type_formdiagram = 'arch'
# discretisation = 20
# folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# csv_file = os.path.join(folder, title + '_data.csv')
# print(csv_file)
# size_parameters, solutions_min, solutions_max = open_csv(csv_file)
# xs.append(size_parameters)
# mins.append(solutions_min)
# maxs.append(solutions_max)
# legends.append('Semicircular arch')

# diagram_file = os.path.join(folder_main, title_bundled + '_diagramA.pdf')
# colors = ['#D2B4DE', '#F39C12', '#2ECC71', '#17202A'][:len(legends)]
# print(diagram_file)
# diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.15, 0.02], [130, 10]]).show()
