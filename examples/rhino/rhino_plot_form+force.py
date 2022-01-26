from compas_tno.diagrams import FormDiagram, ForceDiagram
from compas_tno.rhino import FormArtist, ForceArtist
from compas_tno.rhino import FormObject

from compas_rhino.artists import MeshArtist

import rhinoscriptsyntax as rs

import compas_tno

FORM = []
FORCE = []

FORM.append('/Users/mricardo/compas_dev/compas_tno/data/form-minthk.json')
FORM.append('/Users/mricardo/compas_dev/compas_tno/data/form-minthrust.json')
FORM.append('/Users/mricardo/compas_dev/compas_tno/data/form-maxthrust.json')
FORM.append('/Users/mricardo/compas_dev/compas_tno/data/form-compl.json')

FORCE.append('/Users/mricardo/compas_dev/compas_tno/data/force-minthk.json')
FORCE.append('/Users/mricardo/compas_dev/compas_tno/data/force-minthrust.json')
FORCE.append('/Users/mricardo/compas_dev/compas_tno/data/force-maxthrust.json')
FORCE.append('/Users/mricardo/compas_dev/compas_tno/data/force-compl.json')

for i in [3]:
    form = FormDiagram.from_json(FORM[i])
    force = ForceDiagram.from_json(FORCE[i])
    
    pz = 0.0

    formartist = FormArtist(form)
    print(formartist)
    
    forceartist = ForceArtist(force)
    print(forceartist)
    
    for key in form.vertices():
        pz += form.vertex_attribute(key, 'pz')
        
    print('Type: ', i, pz)

    formartist.pipes_scale = 0.0005
    formartist.layer = 'FormDiagram-' + str(i)

    formartist.draw_forcepipes(tol=1e-4)
    formartist.draw_reactions(scale=0.01)
    formartist.draw_cracks()
    
    #formartist.draw_mesh()
    forceartist.layer = 'ForceDiagram-' + str(i)
    
    forceartist.draw_mesh()
    # used by hand scale of 0.05 in the mesh

    rs.EnableRedraw(True)
