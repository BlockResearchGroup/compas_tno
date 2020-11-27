from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import view_solution
from compas_tno.plotters import plot_form


# Discretisation 10:

# Cross_fd

max_n = '/Users/mricardo/compas_dev/me/max_n/crossvault/cross_fd/crossvault_cross_fd_discr_10_offset-method_min_thk_n_0.25190721006996786.json'
min_t = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_10_min_thk_t_0.3129621556711771.json'

form_n = FormDiagram.from_json(max_n)
form_t = FormDiagram.from_json(min_t)

# Discretisation 14:

# Cross_fd

max_n = '/Users/mricardo/compas_dev/me/max_n/crossvault/cross_fd/crossvault_cross_fd_discr_14_offset-method_min_thk_n_0.2705731208786158.json'
min_t = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_10_min_thk_t_0.3129621556711771.json'

