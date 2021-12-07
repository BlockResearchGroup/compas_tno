
from compas.datastructures import mesh_smooth_centroid
from compas.datastructures import mesh_smooth_area
from compas.datastructures import mesh_smooth_centerofmass

from compas.geometry import closest_point_on_line

from compas_tno.algorithms.equilibrium import equilibrium_fdm

__all__ = [
    'constrained_smoothing',
    'apply_sag'
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

    Return
    ----------
    mesh: Mesh
        The smoothed mesh.
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

    return mesh


def apply_sag(form, boundary_force=10.0, signe_compression=-1.0):  # probably move location
    """Relax the mesh with FDM assuming higher force in the boudary elements.

    Parameters
    ----------
    form : FormDiagram
        The formdiagram to apply sag.
    boundary_force : float
        Force density in the edges on the boundary.
        The default value is ``10.0``.
    signe_compression : float
        Indicate the default sign for compression.
        The default value is ``-1.0``.

    Return
    ----------
    form: FormDiagram
        The relaxed form diagram.

    """

    for u, v in form.edges():
        form.edge_attribute((u, v), 'q', signe_compression*1.0)

    for u, v in form.edges_on_boundary():
        form.edge_attribute((u, v), 'q', signe_compression*boundary_force)

    equilibrium_fdm(form)

    for key in form.vertices():
        form.vertex_attribute(key, 'z', 0.0)

    return form
