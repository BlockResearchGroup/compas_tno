from compas_tno import analysis
from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
import math
import csv
import os

form = FormDiagram.create_cross_form()
vault = Shape.create_pointedcrossvault()

R_over_L = 0.2
lambd = 0.8
thk_0 = 0.5
discr = 14

lambds = [1.0]
R_over_Ls = [0.794]

results = {}

# for lambd in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
for lambd in lambds:
    for spr_angle in [20]:
        for R_over_L in R_over_Ls:

            folder = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/parametric_fd/lambd_{}/spr_{}'.format(lambd, spr_angle)
            # os.makedirs(folder, exist_ok=True)
            title = 'parametric_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}.json'.format(discr, lambd, spr_angle, R_over_L)
            load_form = os.path.join(folder, title)

            img = 'parametric_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}.png'.format(discr, lambd, spr_angle, R_over_L)
            save_img = os.path.join(folder, img)

            form = FormDiagram.from_json(load_form)

            shape_nice = Shape.from_formdiagram_and_attributes(form)

            print(save_img)

            view: Viewer = Viewer(form, show_grid=False)
            view.settings['scale.edge.thk_absolute'] = view.settings['scale.edge.thk_absolute'] * 10
            view.draw_thrust()
            view.draw_shape()
            view.draw_cracks()
            view.show()
