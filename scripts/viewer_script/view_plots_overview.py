from compas_tno.diagrams import FormDiagram

from compas_tno.viewers import Viewer

# from compas_tno.viewers import Viewer
# from compas_tno.plotters import plot_form
# from compas_tno.algorithms import compute_reactions

address = '/Users/mricardo/compas_dev/me/compl_energy/assym/dome/radial_fd/sign_-1/dome_radial_fd_discr_[20, 16]_Ecomp-linear_thk_50.0.json'

form = FormDiagram.from_json(address)

# plot_form(form).show()

view = Viewer(form)
view.show_solution()
