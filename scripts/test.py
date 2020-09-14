
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
import os
import string


# address = '/Users/mricardo/compas_dev/me/shape_comparison/domicalvault/cross_fd/domicalvault_cross_fd_discr_10_min_thk_50.0.json'

type_structure = 'pointed_crossvault'
# type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
type_formdiagram = 'topology-mix'  # Try also 'fan_fd'
discretisation = 10

# hc_list = [6.32, 6.71, 7.42, 7.75, 8.06]
hc = 6.32

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'camb_h='+str(hc))


files = os.listdir(folder)
# print(files)
jsons = []
for f in files:
    extension = f.split('.')[-1]
    if extension == 'json':
        jsons.append('.'.join(f.split('.')[:-1]))

files_dict = {}

for title in jsons:

    title_split = title.split('_')
    thk = float(title_split[-1])
    type_structure = '_'.join(title_split[:2])
    type_formdiagram = title_split[2]

    if type_formdiagram == 'cross' or type_formdiagram == 'fan':
        type_formdiagram = type_formdiagram + '_fd'

    if 'max' in title_split:
        objective = 'max'
    if 'min' in title_split:
        objective = 'min'
    if 'sag' in title:
        sag = int(title.split('sag_')[-1].split('_')[0])
    else:
        sag = False
    if 'smooth' in title:
        smooth = True
    else:
        smooth = False

    data_file = {
        'thk': thk,
        'type_structure': type_structure,
        'type_formdiagram': type_formdiagram,
        'objective': objective,
        'sag': sag,
        'smooth': smooth,
    }

    files_dict[title] = data_file


filters = {'smooth': True}


limit_thk_min = 100
limit_thk_max = 100
limit_form_min = None
limit_form_max = None

for key, values in files_dict.items():
    proceed = True
    if filters:
        for filter in filters:
            if filters[filter] == values[filter]:
                pass
            else:
                proceed = False
                break
    if proceed:
        if values['objective'] == 'min':
            if values['thk'] < limit_thk_min:
                limit_thk_min = values['thk']
                limit_form_min = {key: values}
        if values['objective'] == 'max':
            if values['thk'] < limit_thk_max:
                limit_thk_max = values['thk']
                limit_form_max = {key: values}

print(limit_form_min)
print(limit_form_max)
