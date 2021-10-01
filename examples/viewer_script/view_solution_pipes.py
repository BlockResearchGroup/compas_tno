from compas_tno.diagrams import FormDiagram

from compas_tno.viewers import view_solution2

address = '/Users/mricardo/compas_dev/me/general_opt/min_thk/crossvault/cross_fd/mov_c_0.1/crossvault_cross_fd_discr_10_min_thk_50.0.json'

form = FormDiagram.from_json(address)

view_solution2(form).show()
