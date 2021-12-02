from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
from compas_tno.shapes import Shape
from compas_tno.rhino import ShapeArtist
import rhinoscriptsyntax as rs

import compas
print(compas.__version__)
print(compas)

#for i, j, lambd in [[2, 6, 1], [2, 6, 0.75], [2, 6, 0.50]]:  # [0, 2, 6]
for i, j, lambd in [[2, 6, 0.40]]:  # [0, 2, 6]
    
    address = '/Users/mricardo/compas_dev/compas_tno/data/form_q={}_{}-lambda_{}.json'.format(i, j, lambd)
    
    print(address)
    
    form = FormDiagram.from_json(address)
    
    for edge in form.edges():
        q = form.edge_attribute(edge, 'q')
        form.edge_attribute(edge, 'q', -1 * q)
        
    artist = FormArtist(form)

    artist.pipes_scale = 0.02

    artist.draw_forcepipes(tol=1e-4)
    #artist.draw_reactions(scale=0.005)
    #artist.draw_cracks()
    artist.draw_thrust()

    rs.EnableRedraw(True)