import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

xs = []
mins = []
maxs = []
legends = []

folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'cut_arch')
title_bundled = 'Damage_study'

csv_file = os.path.join(folder_main, 'csv', 'diagram.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
print(size_parameters)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Initial Arch')

diagram_file = os.path.join(folder_main, title_bundled + '_diagramA.pdf')
colors = ['#17202A', '#F39C12', '#2ECC71'][:len(legends)]
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.20, 0.10], [60, 30]]).show()


csv_file = os.path.join(folder_main, 'csv', 'diagram_2.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Cuts_2')

diagram_file = os.path.join(folder_main, title_bundled + '_diagramB.pdf')
colors = ['#17202A', '#F39C12', '#2ECC71'][:len(legends)]
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.20, 0.10], [60, 30]]).show()


csv_file = os.path.join(folder_main, 'csv', 'diagram_3.csv')
print(csv_file)
size_parameters, solutions_min, solutions_max = open_csv(csv_file)
xs.append(size_parameters)
mins.append(solutions_min)
maxs.append(solutions_max)
legends.append('Cuts_1')

diagram_file = os.path.join(folder_main, title_bundled + '_diagramC.pdf')
colors = ['#17202A', '#F39C12', '#2ECC71'][:len(legends)]
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=False, xy_limits=[[0.20, 0.10], [60, 30]]).show()
