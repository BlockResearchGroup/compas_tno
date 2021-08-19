from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

address_shape = '/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=0/crossvault_fan_fd_discr_14_deg=0_min_thk_50.0.json'

address_shape = '/Users/mricardo/compas_dev/me/images/pavilion3.json'


form = FormDiagram.from_json(address_shape)
artist = FormArtist(form)
artist.scale_forces = 1.0

#artist.draw_edges()
#artist.draw_intrados()
#artist.draw_middle()
artist.draw_thrust()
artist.draw_forcepipes()
#artist.draw_extrados()

rs.EnableRedraw(True)
