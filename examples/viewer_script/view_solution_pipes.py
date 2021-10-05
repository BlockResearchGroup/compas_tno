from compas_tno.diagrams import FormDiagram

from compas_tno.viewers import view_solution2
from compas_tno.plotters import plot_form
from compas_tno.algorithms import reactions

# address = '/Users/mricardo/compas_dev/me/general_opt/min_thk/crossvault/cross_fd/mov_c_0.1/crossvault_cross_fd_discr_10_min_thk_50.0.json'
address = '/Users/mricardo/compas_dev/me/compl_energy/assym/crossvault/cross_fd/mov_c_0.1/sign_1/crossvault_cross_fd_discr_10_Ec_thk_50.0.json'

form = FormDiagram.from_json(address)
reactions(form)

plot_form(form).show()
view_solution2(form).show()
