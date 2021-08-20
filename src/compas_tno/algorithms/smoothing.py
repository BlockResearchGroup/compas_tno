
from compas.datastructures import mesh_smooth_centroid
from compas.datastructures import mesh_smooth_area
from compas.datastructures import mesh_smooth_centerofmass

from compas.geometry import closest_point_on_line

__all__ = [
    'constrained_smoothing'
]


def constrained_smoothing(mesh, kmax=100, damping=0.5,  constraints={}, algorithm='centroid'):
    """Constrained smoothing of a mesh. Constraints can be points and lines.

    Parameters
    ----------
    mesh : Mesh
        A mesh to smooth.
    kmax : int
        Number of iterations for smoothing. Default value 100.
    damping : float
        Damping value for smoothing between 0 and 1. Default value 0.5.
    constraints : dict
        Dictionary of constraints as vertex keys pointing to fixed points or lines.
    algorithm : string
        Type of smoothing algorithm to apply (classic centroid or area-based). Classic centroid by default.
    """

    def callback(k, args):

        mesh, constraints = args

        for vkey, constraint in constraints.items():
            if constraint is None:
                continue
            elif len(constraint) == 3:  # It is a Point
                x, y, z = constraint
            elif len(constraint) == 2:  # It is a Line
                x, y, z = closest_point_on_line(mesh.vertex_coordinates(vkey), constraint)
            else:
                continue

            mesh.vertex[vkey]['x'] = x
            mesh.vertex[vkey]['y'] = y
            mesh.vertex[vkey]['z'] = z

    func = {'centroid': mesh_smooth_centroid, 'area': mesh_smooth_area, 'centerofmass': mesh_smooth_centerofmass}

    if algorithm not in func:
        algorithm = 'centroid'

    func[algorithm](mesh, kmax=kmax, damping=damping, callback=callback, callback_args=[mesh, constraints])
