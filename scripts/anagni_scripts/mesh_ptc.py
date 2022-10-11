# Script to mesh the pointcloud of the intrados of the S. Agostino
# dome, a case study of the Anagni Summer school.

import open3d as o3d
import time

# read point cloud from files
start_time = time.time()
pcd_file = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/AgostinoDEF.pts'
# pcd_file = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_downsized_25.0.pts'

pcd = o3d.io.read_point_cloud(pcd_file)

# select voxel size to downsample the pointcloud
voxel_size = 0.10
downpcd = pcd.voxel_down_sample(voxel_size=voxel_size)

# print number of vertices
print('Original data', pcd)
print('Voxel down sample', downpcd)


pcd_downsized = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_downsized_' + str(voxel_size*100) + '.pts'
mesh_downsized = '/Users/mricardo/Documents/ETH/Teaching/2021_Anagni/SAgostino/Agostino_mesh_' + str(voxel_size*100) + '.obj'

# o3d.io.write_point_cloud(pcd_downsized, downpcd)
# print('Written the downsized in:', pcd_downsized)

end_time = time.time() - start_time

print('Time elapsed to downsized', end_time)
print("Recompute the normal of the downsampled point cloud")

downpcd.estimate_normals(
    search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))

downpcd.orient_normals_consistent_tangent_plane(100)
o3d.visualization.draw_geometries([downpcd])

o3d.visualization.draw_geometries([downpcd])

with o3d.utility.VerbosityContextManager(
        o3d.utility.VerbosityLevel.Debug) as cm:
    mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
        downpcd, depth=9)

print(mesh)

o3d.io.write_triangle_mesh(mesh_downsized, mesh)

o3d.visualization.draw_geometries([mesh])
