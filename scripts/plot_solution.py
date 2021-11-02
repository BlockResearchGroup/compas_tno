from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram

type_formdiagram = 'cross_fd'
span = 10.0
discretisation = 14
json = '/Users/mricardo/compas_dev/me/compl_energy/crossvault/cross_fd/mov_c_0.1/corner/sign_1/crossvault_cross_fd_discr_14_Ecomp-linear_thk_50.0.json'

form = FormDiagram.from_json(json)

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form_base = FormDiagram.from_library(data_diagram)

plot_superimposed_diagrams(form, form_base).show()

