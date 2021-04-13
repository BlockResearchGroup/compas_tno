from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

address_shape = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/tributary_area/R=5.0/deg=0/primal_pointed_crossvault_topology-crossbraced_discr_14.json'

shapes = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.1147/min_thk/deg=20/pointed_crossvault_cross_fd_discr_20_shape_thk_50.0.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=7.08/min_thk/deg=20/pointed_crossvault_cross_fd_discr_20_shape_thk_50.0.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=7.9422/min_thk/deg=20/pointed_crossvault_cross_fd_discr_20_shape_thk_50.0.json']

address_shape = shapes[2]

form = FormDiagram.from_json(address_shape)
artist = FormArtist(form)

#artist.draw_edges()
artist.draw_intrados()
#artist.draw_middle()
artist.draw_extrados()

rs.EnableRedraw(True)
