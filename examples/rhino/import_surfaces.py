from compas_rhino.conversions import RhinoSurface
from compas_rhino.conversions import RhinoMesh
from compas_rhino.utilities import select_surface
from compas_rhino.utilities import select_mesh
import compas_rhino

compas_rhino.rs.UnselectAllObjects()

jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_rectangular/'

guid_intra = select_surface('Select the Intrados')
# surface = RhinoSurface.from_guid(guid)

# compas_surf = surface.to_compas()
# compas_surf = surface.to_compas_mesh()
# compas_surf = surface.to_compas_quadmesh()

# compas_surf.to_json(jsonpath + '1-intrados-mesh.json')

compas_rhino.rs.UnselectAllObjects()

guid_extra = select_surface('Select the Extrados')
# surface = RhinoSurface.from_guid(guid)

# compas_surf = surface.to_compas()
#compas_surf = surface.to_compas_mesh()
# compas_surf = surface.to_compas_quadmesh()

# compas_surf.to_json(jsonpath + '1-extrados-mesh.json')

compas_rhino.rs.UnselectAllObjects()

guid_mesh = select_mesh('Select mesh base')
mesh = RhinoMesh.from_guid(guid_mesh).to_compas()

print(mesh)

compas_rhino.rs.UnselectAllObjects()

points =  mesh.vertices_attributes('xyz')
    
points_intra = compas_rhino.rs.ProjectPointToSurface(points, guid_intra, (0,0,1))
points_extra = compas_rhino.rs.ProjectPointToSurface(points, guid_extra, (0,0,1))

# xxxx
for i, vertex in enumerate(mesh.vertices()):
    mesh.vertex_attribute(vertex, 'z', points_intra[i][2])

path = jsonpath + '1-intrados-mesh.json'
mesh.to_json(path)
print('mesh saved to:', path)

# xxxx
for i, vertex in enumerate(mesh.vertices()):
    mesh.vertex_attribute(vertex, 'z', points_extra[i][2])

path = jsonpath + '1-extrados-mesh.json'
mesh.to_json(path)
print('mesh saved to:', path)