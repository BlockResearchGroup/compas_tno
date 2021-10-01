from compas.datastructures import Mesh
from compas.geometry import Pointcloud
from compas_plotters import MeshPlotter

mesh_file = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_mesh_downsized.ply'

pc = Pointcloud.from_ply(mesh_file)

mesh = Mesh.from_ply(mesh_file)

plotter = MeshPlotter(mesh)
plotter.draw_edges()
plotter.draw_vertices()
plotter.show()
