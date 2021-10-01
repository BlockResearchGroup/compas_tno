from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

address_shape = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-fansmooth/R=11.5/min_thk/pointed_crossvault_topology-fansmooth_discr_14_lp_thk_50.0.json'

address_shape = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/tributary_area/R=5.0/deg=0/dual_pointed_crossvault_cross_fd_discr_14.json'

address_shape = '/Users/mricardo/compas_dev/me/minmax/dome/radial_spaced/radial_spaced_discr_8_20_min_t=207.json'


form = FormDiagram.from_json(address_shape)
artist = FormArtist(form)

artist.draw()

#artist.draw_edges()
#artist.draw_intrados()
#artist.draw_middle()
#artist.draw_extrados()
#artist.draw_extrados()

rs.EnableRedraw(True)
