import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

discretisation = 10
type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'
xs = []
mins = []
maxs = []
legends = []

folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)

for ratio in [0.10, 0.05]:

    title_ratio = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_partial_ratio_' + str(ratio)
    csv_file = os.path.join(folder, title_ratio + '_data.csv')
    print(csv_file)
    size_parameters, solutions_min, solutions_max = open_csv(csv_file)
    # diagram_of_thrust(size_parameters, solutions_min, solutions_max, limit_state=True).show()
    xs.append(size_parameters)
    mins.append(solutions_min)
    maxs.append(solutions_max)
    legends.append('Hor. React. ' + str(100*ratio)+ '%')


for fill_percentage in [1.0, 0.75]:

    title_fill = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_fill_' + str(fill_percentage)
    csv_file = os.path.join(folder, title_fill + '_data.csv')
    print(csv_file)
    size_parameters, solutions_min, solutions_max = open_csv(csv_file)
    # diagram_of_thrust(size_parameters, solutions_min, solutions_max, limit_state=True).show()
    xs.append(size_parameters)
    mins.append(solutions_min)
    maxs.append(solutions_max)
    legends.append('Infill ' + str(100*fill_percentage)+ '%')

title_normal = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
csv_file = os.path.join(folder, title_normal + '_data.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
# diagram_of_thrust(size_parameters, solutions_min, solutions_max, limit_state=True).show()
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Initial')

title_bundled = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_all_'
diagram_file = os.path.join(folder, title_bundled + '_diagram.pdf')
print(diagram_file)
colors = ['#D2B4DE', '#D2B4DE', '#F39C12', '#F39C12', '#17202A']
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, xy_limits=[[0.10, 0.035], [70, 30]]).show()
