from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.rhino import FormArtist
from compas_tno.rhino import ForceArtist
import compas_rhino

plot_form = False
plot_force = True

path = '/Users/mricardo/compas_dev/compas_tno/data/form.json'

path = '/Users/mricardo/compas_dev/me/compl_energy/dome/outwards/radial_fd/dome_radial_fd_discr_[16, 20]_Ecomp-linear_thk_50.0.json'
#path = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/radial_fd/dome_radial_fd_discr_[16, 20]_Ecomp-linear_thk_50.0.json'
path = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/parametric/form_sigularity_ni_2_ns_5.json/circular_sigularity_ni_2_ns_5.json_Ecomp-linear_thk_50.0.json'
#path = '/Users/mricardo/compas_dev/me/compl_energy/crossvault/corner/cross_fd/crossvault_cross_fd_discr_14_Ecomp-linear_thk_50.0.json'
#path = '/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open/cross_fd/pavillionvault_cross_fd_discr_14_spr_30.0_Ecomp-linear_thk_50.0.json'
#path = '/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open/cross_fd/hor_loads/pavillionvault_cross_fd_discr_14_spr_30.0_Ecomp-linear_thk_50.0.json'

ni = 4
ns = 2
code = 'ni_{}_ns_{}'.format(ni, ns)
path = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/parametric/form_sigularity_' + code + '.json/circular_sigularity_' + code + '.json_Ecomp-linear_thk_50.0.json'.format(ni, ns)

if plot_form:

    form = FormDiagram.from_json(path)

    artist = FormArtist(form)
    artist.draw_thrust()
    artist.draw_cracks()
    #artist.draw_reactions(scale=0.005)  # crossvault
    artist.draw_reactions(scale=0.05)  # dome
    artist.redraw()

# path = '/Users/mricardo/compas_dev/me/compl_energy/crossvault/corner/cross_fd/crossvault_cross_fd_discr_14_Ecomp-linear_thk_50_force.json'
# path = '/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open/cross_fd/pavillionvault_cross_fd_discr_14_spr_30_force.json'

path = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/parametric/form_sigularity_' + code + '_force.json'

force = ForceDiagram.from_json(path)

#from compas.datastructures import Mesh
#from compas.geometry import Scale
#s = Scale.from_factors([0.1, 0.1, 0.1])
#force.transform(s)

if plot_force:
    
    artist = ForceArtist(force)
    artist.draw()

compas_rhino.redraw()