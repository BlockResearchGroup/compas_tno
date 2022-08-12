
from jax.numpy import dot
from jax.numpy import hstack
from jax.numpy.linalg import norm
# from jax import grad

from jax.numpy import cross
from jax.numpy import transpose


def weights_from_xyz_jax(xyz, F, V0, V1, V2, thk=0.5, density=20.0):

    v0 = dot(V0, xyz)
    v1 = dot(V1, xyz) - v0
    v2 = dot(V2, dot(F, xyz)) - v0
    ag = norm(cross(v1, v2), axis=1)

    pz = 0.25 * thk * density * dot(transpose(V0), ag)

    return pz


def weights_from_xyz_separate_jax(x, y, z, F, V0, V1, V2, thk=0.5, density=20.0):

    xyz = hstack([x, y, z])

    v0 = dot(V0, xyz)
    v1 = dot(V1, xyz) - v0
    v2 = dot(V2, dot(F, xyz)) - v0
    ag = norm(0.25 * cross(v1, v2), axis=1)

    pz = 0.25 * thk * density * dot(transpose(V0), ag)

    return pz
