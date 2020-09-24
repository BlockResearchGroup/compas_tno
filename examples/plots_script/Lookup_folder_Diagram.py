import os
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import save_csv_row
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import lookup_folder
from compas_tno.plotters import filter_min_thk
from compas_tno.plotters import plot_form
from compas_tno.diagrams import FormDiagram

type_structure = 'dome'
type_formdiagram = 'radial_spaced_fd'  # y also 'fan_fd'
discretisation = [8, 20]
hc = None


if hc:
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
else:
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# ----------------------- Lookup Folder --------------------------

files_dict = lookup_folder(folder)
limit_form_min, limit_form_max = filter_min_thk(files_dict, filters={'smooth': True})
print(limit_form_min)
print(limit_form_max)

sol_min = []
sol_max = []
thk_min = []
thk_max = []

# If need to filter see the dparent function
for key, values in files_dict.items():
    if values['objective'] == 'min':  # and values['thk']<=40:
        thk_min.append(values['thk']/100)
        form_address = os.path.join(folder, key + '.json')
        form = FormDiagram.from_json(form_address)
        sol_min.append(form.thrust()/form.lumped_swt())
        if values['thk'] == 50.0:
            print('min')
            plot_form(form, show_q=False, cracks=True).show()
    if values['objective'] == 'max':  # and values['thk']<=40:
        thk_max.append(values['thk']/100)
        form_address = os.path.join(folder, key + '.json')
        form = FormDiagram.from_json(form_address)
        sol_max.append(-1*form.thrust()/form.lumped_swt())
        if values['thk'] == 50.0:
            print('max')
            plot_form(form, show_q=False, cracks=True).show()

zipped_lists = zip(thk_min, sol_min)
sorted_zipped_lists = sorted(zipped_lists, reverse=True)
thk_min = [element for element, _ in sorted_zipped_lists]
sol_min = [element for _, element in sorted_zipped_lists]

zipped_lists = zip(thk_max, sol_max)
sorted_zipped_lists = sorted(zipped_lists, reverse=True)
thk_max = [element for element, _ in sorted_zipped_lists]
sol_max = [element for _, element in sorted_zipped_lists]

# print(sorted_list1)
# print(sorted_list2)

# print(files_dict)
# print('-'*10)
# print(sol_min)
# print(thk_min)
# print('-'*10)
# print(sol_max)
# print(thk_max)

thicknesses = [thk_min, thk_max]
solutions = [sol_min, sol_max]
print('\nThicknesses:')
print('-'*10)
print(thicknesses)
print('\nSolutions:')
print('-'*10)
print(solutions)

# ----------------------- Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
# thicknesses, solutions = open_csv_row(csv_file)

# resave with different name
save_csv_row(thicknesses, solutions, path=csv_file, title=type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation))

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True).show()

xy_limits = [[0.50, 0.15], [100, 0]]
img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits).show()
