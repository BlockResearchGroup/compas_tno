from compas_tno.diagrams import FormDiagram
from compas_tno.optimisers import Optimiser
from compas_tno.algorithms import reactions
import os
import compas_tno
import math

discretisation = [4, 12]
load_mult = 0.00
load_increase = 0.01
direction_loads = 'px'
optimiser = Optimiser()
optimiser.settings['objective'] = 'max'
type_formdiagram = 'radial_fd'
style_diagonals = 'straight'
# title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + type_formdiagram + '_' + style_diagonals + '_' + optimiser.settings['objective']
folder = compas_tno.get('/dome/') # Folder to Save the structure
size_parameters = []
size_min = []

while load_mult <= 0.30:
    title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + type_formdiagram + '_' + style_diagonals + '_' + optimiser.settings['objective']
    # title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective']
    jsonpath = os.path.join(folder, title + '.json')
    form = FormDiagram.from_json(jsonpath)
    reactions(form, plot=False)
    # print(jsonpath)
    thrust = 0
    swt = 0
    for key in form.vertices_where({'is_fixed': True}):
        rx = form.vertex_attribute(key, '_rx')
        ry = form.vertex_attribute(key, '_ry')
        r = math.sqrt(rx**2 + ry**2)
        thrust += r
    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        swt += pz
    # print(thrust)
    # print(swt)
    # print(thrust/swt)
    print('load mult and T/W:', load_mult, thrust, swt, thrust/swt)
    size_parameters.append(load_mult)
    size_min.append(thrust/swt)
    load_mult = round(load_mult + load_increase, 4)

print(size_parameters)
print(size_min)
