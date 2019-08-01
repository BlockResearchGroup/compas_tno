from __future__ import print_function
from __future__ import absolute_import
from __future__ import division


from compas.utilities import pairwise
from compas.geometry import centroid_points
from compas.geometry import add_vectors, scale_vector, length_vector
from compas_tessellation.utilities import i_to_blue, i_to_red, i_to_orange


__all__ = ['AssemblyView', 'BlockView', 'InterfaceView']


class AssemblyView(object):

    def __init__(self, assembly):
        self._assembly = None
        self._xyz = None
        self._vertices = None
        self._is_support = None
        self.assembly = assembly

    # TODO: implement assembly view model, error here. never reached
    @property
    def xyz(self):
        return self._xyz

    @property
    def vertices(self):
        return self.assembly.vertices()

    @property
    def is_support(self):
        return self._is_support

    @property
    def edges(self):
        key_index = self.assembly.key_index()
        for u, v in self.assembly.edges():
            yield key_index[u], key_index[v]

    @property
    def assembly(self):
        return self._assembly

    @assembly.setter
    def assembly(self, assembly):
        self._assembly = assembly
        self._is_support = assembly.get_vertices_attribute('is_support')
        self._xyz = assembly.get_vertices_attributes('xyz')


class BlockView(object):
    
    def __init__(self, block):
        self._block = None
        self._xyz = None
        self._vertices = None
        self._faces = None
        self.block = block

    @property
    def xyz(self):
        return self._xyz

    @property
    def vertices(self):
        return self.block.vertices()

    @property
    def faces(self):
        return self._faces

    @property
    def edges(self):
        key_index = self.block.key_index()
        for u, v in self.block.edges():
            yield key_index[u], key_index[v]

    @property
    def block(self):
        return self._block

    @block.setter
    def block(self, block):
        self._block = block

        key_index = block.key_index()

        xyz = block.get_vertices_attributes('xyz')
        faces = []
        # print(xyz)
        for fkey in block.faces():
            fvertices = [key_index[key] for key in block.face_vertices(fkey)]

            f = len(fvertices)
            if f < 3:
                pass
            elif f == 3:
                faces.append(fvertices)
            elif f == 4:
                a, b, c, d = fvertices
                faces.append([a, b, c])
                faces.append([c, d, a])
            else:
                o = block.face_centroid(fkey)
                v = len(xyz)
                xyz.append(o)
                for a, b in pairwise(fvertices + fvertices[0:1]):
                    faces.append([a, b, v])

        self._xyz = xyz
        self._faces = faces


class InterfaceView(object):

    def __init__(self, interface):
        self._interface = None
        self._xyz = None
        self._faces = None
        self._colors = None
        self._frictions = None
        self._forces = None
        self.interface = interface

    @property
    def xyz(self):
        return self._xyz

    @property
    def faces(self):
        return self._faces

    @property
    def colors(self):
        return self._colors

    @property
    def frictions(self):
        return self._frictions

    @property
    def forces(self):
        return self._forces

    @property
    def interface(self):
        return self._interface

    @interface.setter
    def interface(self, interface):
        self._interface = interface

        faces = []
        xyz = interface['interface_points']
        w = interface['interface_uvw'][2]
        fs = interface['interface_forces']
        scale = 1
        eps = 1e-2

        if fs is None or w is None:
            print("error")
            return

        v = len(xyz)

        if v < 3:
            pass

        elif v == 3:
            faces.append([0, 1, 2])

        elif v == 4:
            faces.append([0, 1, 2])
            faces.append([2, 3, 0])
        else:
            c = centroid_points(xyz)
            xyz.append(c)
            for a, b in pairwise(list(range(0, v))):
                faces.append([a, b, v])
            faces.append([b, 0, v])

        forces = [f['c_np'] - f['c_nn'] for f in fs]
        fmin = min(forces)
        fmax = max(forces)

        # compute friction color
        u, v = interface['interface_uvw'][0], interface['interface_uvw'][1]
        fforces = [length_vector(add_vectors(scale_vector(u, f['c_u']), scale_vector(v, f['c_v']))) for f in fs]
        ffmin = min(fforces)
        ffmax = max(fforces)

        # coloring interface
        rgb = []
        lines = []
        friction_rgb = []

        def blue(i):
            i = max(i, 0.0)
            i = min(i, 1.0)
            rb = min((1 - i) * 255, 255)
            return (255, 255, int(rb))

        for i in range(len(xyz)):
            if i >= len(fs):
                continue

            # print('number of vertices', i)
            sp = xyz[i]
            f = forces[i]
            if f > eps:
                force = 'compression'
                color = (0.0, 0.0, 1.0)
                rgb.append(i_to_blue(abs(1 * f / fmax)))
                # rgb.append((1-(abs(1 * f / fmax)), 1-(abs(1 * f / fmax)), 1))
            elif f < -eps:
                force = 'tension'
                color = (1.0, 0.0, 0.0)
                rgb.append(i_to_red(abs(1 * f / fmin)))
                # rgb.append((1, 1-abs(1 * f / fmin), 1-abs(1 * f / fmin)))
            else:
                force = 'zero'
                color = (1.0, 1.0, 1.0)
                rgb.append(color)

            if fforces[i] > eps:
                friction_rgb.append(i_to_orange(abs(1 * fforces[i] / ffmax)))
            else:
                friction_rgb.append((1.0, 1.0, 1.0))

            lines.append({
                # 'start': [sp[axis] - 0.5 * scale * f * w[axis] for axis in range(3)],
                'start': sp,
                'end': [sp[axis] + 1 * scale * f * w[axis] for axis in range(3)],
                'indices': (i * 2, i * 2 + 1),
                'color': color,
                'force': force
            })

        # assign to face center color as average
        if len(xyz) > 4:
            fa = sum(forces) / (len(xyz) - 1)
            ffa = sum(fforces) / (len(xyz) - 1)

            if scale * fa > eps:
                rgb.append(i_to_blue(abs(fa / fmax)))
                # rgb.append((0.0, 0.0, abs(1 * fa / fmax)))
            elif scale * fa < -eps:
                rgb.append(i_to_red(abs(fa / fmin)))
                # rgb.append((abs(1 * fa / fmin), 0.0, 0.0))
            else:
                rgb.append((255, 255, 255))

            if ffa > eps:
                friction_rgb.append(i_to_orange(abs(ffa / ffmax)))
            else:
                friction_rgb.append((1.0, 1.0, 1.0))

        self._forces = lines
        self._colors = rgb
        self._frictions = friction_rgb
        self._xyz = xyz
        self._faces = faces



# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    pass
