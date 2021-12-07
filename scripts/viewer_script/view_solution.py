from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.plotters import plot_form
import os
os.environ['QT_MAC_WANTS_LAYER'] = '1'

# # Discretisation 10:

# # Cross_fd

# max_n = '/Users/mricardo/compas_dev/me/max_n/crossvault/cross_fd/crossvault_cross_fd_discr_10_offset-method_min_thk_n_0.25190721006996786.json'
# min_t = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_10_min_thk_t_0.3129621556711771.json'

# form_n = FormDiagram.from_json(max_n)
# form_t = FormDiagram.from_json(min_t)

# # Discretisation 14:

# # Cross_fd

# max_n = '/Users/mricardo/compas_dev/me/max_n/crossvault/cross_fd/crossvault_cross_fd_discr_14_offset-method_min_thk_n_0.2705731208786158.json'
# min_t = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_10_min_thk_t_0.3129621556711771.json'


# address = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/refined/min_thk/deg=40/pointed_crossvault_cross_fd_discr_14_R=5.326621016707531_min_thk_4.984461856922818.json'

# form = FormDiagram.from_json(address)
# view = Viewer(form)
view.show_solution()

# address = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/refined/min_thk/deg=40/pointed_crossvault_fan_fd_discr_14_R=6.239890208210451_min_thk_11.196091800459316.json'
address = '/Users/mricardo/compas_dev/me/min_thk/dome/radial_fd/min_max/dome_radial_fd_discr_[20, 16]_max_thk_49.99999999999999.json'

form = FormDiagram.from_json(address)

pzt = 0.0
for key in form.vertices():
    pz = form.vertex_attribute(key, 'pz')
    pzt += pz

fmax = 0.0
qmax = 0.0
for u, v in form.edges():
    q = form.edge_attribute((u, v), 'q')
    l = form.edge_length(u, v)
    f = q * l
    if q > qmax:
        qmax = q
    if f > fmax:
        fmax = f

scalefac = 1570/pzt
print('Maximum q:', qmax)
print('Maximum f:', fmax)
print('Total pz:', pzt)
print('Total pz*20:', pzt*20)
print('Scale Factor:', scalefac)
print('Maximum q and f scaled:', scalefac*qmax, scalefac*fmax)


# view = Viewer(form)
view.show_solution()
