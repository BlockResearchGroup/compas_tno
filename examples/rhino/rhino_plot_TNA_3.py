from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
from compas_tno.shapes import Shape
from compas_tno.rhino import ShapeArtist
import rhinoscriptsyntax as rs


address = '/Users/mricardo/compas_dev/compas_tno/data/form_tna_cap.json'

form = FormDiagram.from_json(address)
        
artist = FormArtist(form)

artist.pipes_scale = 0.0001

#artist.draw_forcepipes(tol=1e-4)
artist.draw_reactions(scale=0.005)
#artist.draw_cracks()
#artist.draw_thrust()

rs.EnableRedraw(True)