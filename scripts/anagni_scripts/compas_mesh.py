from compas.datastructures import Mesh
from compas.geometry import Pointcloud
# from compas_plotters import MeshPlotter
from compas_view2 import app

mesh_file = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_downsized_10.0.pts'

pc = Pointcloud.from_ply(mesh_file)

# mesh = Mesh.from_ply(mesh_file)

viewer = app.App()
viewer.add(pc)
viewer.show()

# plotter = MeshPlotter(mesh)
# plotter.draw_edges()
# plotter.draw_vertices()
# plotter.show()
