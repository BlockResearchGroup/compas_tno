import compas_tno
import os
import copy
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import plot_sym_inds

# ------------------------------------------------------------------------------------
# ------------------------------ MAKE THE FROM DIAGRAM -------------------------------
# ------------------------------------------------------------------------------------

span = 10.0  # square span for analysis
thk = 0.50

# Basic parameters

type_structure = 'crossvault'
type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
discretisation = 10
# discretisation = [8, 16]

# ----------------------- 1. Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, span]],
    'discretisation': discretisation,
    'fix': 'corners',
    'center': [span/2, span/2, 0.0],
    'radius': span/2,
    'r_oculus': 0.0,
    'diagonal': None,
    'partial_diagonal': None,
}

form = FormDiagram.from_library(data_diagram)

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)

# --------------------- 3. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'MMA'
optimiser.data['solver'] = 'MMA'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = False
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 10e+10


# --------------------- 4. Shape with initial THK ---------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, span]],
    't': 0.0,
    'center': [span/2, span/2, 0.0],
    'radius': span/2,
    'r_oculus': 0.0,
}

vault = Shape.from_library(data_shape)

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_envelope()
analysis.apply_selfweight()
# analysis.apply_symmetry(center_point=[span/2, span/2, 0.0])
analysis.set_up_optimiser()
analysis.run()
# plot_independents(form).show()

plot_form(form, show_q=False, cracks=True).show()

json_path = compas_tno.get('test.4.json')
form.to_json(json_path)


# ------------------------------------------------------------------------------------
# -------------------------- EXTRACT TO SAVE TIME ------------------------------------
# ------------------------------------------------------------------------------------


# json_path = compas_tno.get('test.json')
# form = FormDiagram.from_json(json_path)

# form.apply_symmetry(center=[5.0, 5.0, 0.0])
# Asym = form.build_symmetry_matrix(printout=True)
# plot_sym_inds(form).show()
# import sys
# import numpy
# numpy.set_printoptions(threshold=sys.maxsize)
# print(Asym)


# import matplotlib.pyplot as plt

# # Display matrix
# plt.matshow(Asym)
# plt.show()

# plot_sym_inds(form).show()
# plot_independents(form, show_symmetry=True).show()
# form.build_symmetry_matrix(printout=True)
# form.build_symmetry_matrix_supports(printout=True)
# form.assemble_symmetry_matrix(printout=True)

