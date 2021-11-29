from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs

file = '/Users/mricardo/compas_dev/me/compl_energy/assym/dome/radial_fd/sign_-1/dome_radial_fd_discr_[20, 16]_Ecomp-linear_thk_50.0.json'

form = FormDiagram.from_json(file)
artist = FormArtist(form)

artist.pipes_scale = 0.001

#artist.draw_forcepipes(tol=1e-4)
#artist.draw_reactions(draw_as_pipes=True)
#artist.draw_cracks()
artist.draw_thrust()

rs.EnableRedraw(True)