# Script to mesh the pointcloud of the intrados of the S. Agostino
# dome, a case study of the Anagni Summer school.

import open3d as o3d
import time

# read point cloud from files
# pcd_file = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/AgostinoDEF.pts'
# pcd_downsized = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_downsized.pts'
mesh_downsized = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_mesh_downsized.ply'
mesh_obj = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_mesh_downsized.obj'

mesh = o3d.io.read_triangle_mesh(mesh_downsized)

# print number of vertices
print('Loaded data', mesh)

# o3d.visualization.draw_geometries([mesh])

o3d.io.write_triangle_mesh(mesh_obj, mesh)
