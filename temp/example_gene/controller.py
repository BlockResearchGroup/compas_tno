from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from functools import partial

try:
    import PySide2
except ImportError:
    from PySide import QtCore
    from PySide import QtGui
    import PySide.QtGui as QtWidgets
else:
    from PySide2 import QtCore
    from PySide2 import QtGui
    from PySide2 import QtWidgets

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import compas

from compas.datastructures import Mesh

from compas.geometry import centroid_points

from compas.utilities import hex_to_rgb
from compas.utilities import flatten

from compas_tessellation.geometries.quads import create_assembly
from compas_tessellation.app import core
from .model import AssemblyView, BlockView, InterfaceView


__all__ = ['Controller']


hex_to_rgb = partial(hex_to_rgb, normalize=True)


def flist(items):
    return list(flatten(items))


class Controller(core.controller.Controller):
    settings = core.controller.Controller.settings or {}

    settings['mu:value'] = 0.4
    settings['mu:minval'] = 10.0
    settings['mu:maxval'] = 90.0
    settings['mu:step'] = 1.0
    settings['mu:scale'] = 0.01

    settings['cm:value'] = 0.
    settings['cm:minval'] = -30.0
    settings['cm:maxval'] = 30.0
    settings['cm:step'] = 1.0
    settings['cm:scale'] = 0.01

    settings['delta:value'] = 0.0
    settings['delta:minval'] = 0.0
    settings['delta:maxval'] = 90.0
    settings['delta:step'] = 1.0
    settings['delta:scale'] = 0.01

    settings['edges.width:value'] = 1.0
    settings['edges.width:minval'] = 1
    settings['edges.width:maxval'] = 100
    settings['edges.width:step'] = 1
    settings['edges.width:scale'] = 1.0

    settings['vertices.size:value'] = 5
    settings['vertices.size:minval'] = 1
    settings['vertices.size:maxval'] = 100
    settings['vertices.size:step'] = 1
    settings['vertices.size:scale'] = 1.0

    settings['edges.width:value'] = 1.0
    settings['edges.width:minval'] = 1
    settings['edges.width:maxval'] = 100
    settings['edges.width:step'] = 1
    settings['edges.width:scale'] = 1.0

    settings['vertices.color'] = '#0092d2'
    settings['edges.color'] = '#666666'
    settings['faces.color:front'] = '#cccccc'
    settings['faces.color:back'] = '#ff5e99'

    settings['interfaces.color:front'] = '#cccccc'
    settings['interfaces.color:back'] = '#ff5e99'

    settings['vertices.on'] = True
    settings['edges.on'] = True
    settings['faces.on'] = False
    settings['interfaces.on'] = True
    settings['friction.on'] = False
    settings['massopt.on'] = False

    # settings['vertices.labels.on'] = False
    # settings['edges.labels.on'] = False
    # settings['faces.labels.on'] = False

    settings['camera.elevation:value'] = -10
    settings['camera.elevation:minval'] = -180
    settings['camera.elevation:maxval'] = 0
    settings['camera.elevation:step'] = +1
    settings['camera.elevation:scale'] = +1

    settings['camera.azimuth:value'] = +30
    settings['camera.azimuth:minval'] = -180
    settings['camera.azimuth:maxval'] = +180
    settings['camera.azimuth:step'] = +1
    settings['camera.azimuth:scale'] = +1

    settings['camera.distance:value'] = +10
    settings['camera.distance:minval'] = 0
    settings['camera.distance:maxval'] = +100
    settings['camera.distance:step'] = +1
    settings['camera.distance:scale'] = +1
    settings['camera.distance:delta'] = +0.05

    settings['camera.rotation:delta'] = +0.5
    settings['camera.fov:value'] = 50
    settings['camera.near:value'] = 0.1
    settings['camera.far:value'] = 1000

    def __init__(self, app):
        super(Controller, self).__init__(app)
        self._assembly = None
        self._blocks = None
        self._interfaces = None

    @property
    def view(self):
        return self.app.view

    @property
    def assembly(self):
        return self._assembly

    @assembly.setter
    def assembly(self, assembly):
        # self._assembly = assembly
        self._assembly = AssemblyView(assembly)

        self._blocks = []
        for key, attr in assembly.vertices(True):
            self._blocks.append(BlockView(assembly.blocks[key]))

        self._interfaces = []
        for u, v, attr in assembly.edges(True):
            self._interfaces.append(InterfaceView(attr))

    @property
    def blocks(self):
        return self._blocks

    @property
    def interfaces(self):
        return self._interfaces

    def center_assembly(self):
        xyz = self.assembly.get_vertices_attributes('xyz')
        cx, cy, cz = centroid_points(xyz)
        for key, attr in self.assembly.vertices(True):
            attr['x'] -= cx
            attr['y'] -= cy
            attr['z'] -= cz

    # ==========================================================================
    # constructors
    # ==========================================================================

    # ==========================================================================
    # view
    # ==========================================================================

    # ==========================================================================
    # appearance
    # ==========================================================================
    def slide_mu(self, value):
        self.settings['mu:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def edit_mu(self, value):
        self.settings['mu:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def slide_cm(self, value):
        self.settings['cm:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def edit_cm(self, value):
        self.settings['cm:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def slide_delta(self, value):
        self.settings['delta:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def edit_delta(self, value):
        self.settings['delta:value'] = value
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'], self.settings['massopt.on'])

    def slide_size_vertices(self, value):
        self.settings['vertices.size:value'] = value
        self.view.updateGL()

    def edit_size_vertices(self, value):
        self.settings['vertices.size:value'] = value
        self.view.updateGL()

    def slide_width_edges(self, value):
        self.settings['edges.width:value'] = value
        self.view.updateGL()

    def edit_width_edges(self, value):
        self.settings['edges.width:value'] = value
        self.view.updateGL()

    # ==========================================================================
    # visibility
    # ==========================================================================

    def toggle_interfaces(self, state):
        self.settings['interfaces.on'] = state == QtCore.Qt.Checked
        self.view.updateGL()

    def toggle_faces(self, state):
        self.settings['faces.on'] = state == QtCore.Qt.Checked
        self.view.updateGL()

    def toggle_edges(self, state):
        self.settings['edges.on'] = state == QtCore.Qt.Checked
        self.view.updateGL()

    def toggle_vertices(self, state):
        self.settings['vertices.on'] = state == QtCore.Qt.Checked
        self.view.updateGL()

    def toggle_friction(self, state):
        self.settings['friction.on'] = state == QtCore.Qt.Checked
        self.view.glInit()
        self.view.make_buffers()
        self.view.updateGL()

    def toggle_massopt(self, state):
        self.settings['massopt.on'] = state == QtCore.Qt.Checked
        self.create_floor_assembly(self.settings['mu:value'], self.settings['cm:value'], self.settings['delta:value'],
                                   self.settings['massopt.on'])
        self.view.updateGL()

    # ==========================================================================
    # color
    # ==========================================================================

    def change_vertices_color(self, color):
        self.settings['vertices.color'] = color
        self.view.update_vertex_buffer('vertices.color', self.view.array_vertices_color)
        self.view.updateGL()
        self.app.main.activateWindow()

    def change_edges_color(self, color):
        self.settings['edges.color'] = color
        self.view.update_vertex_buffer('edges.color', self.view.array_edges_color)
        self.view.updateGL()
        self.app.main.activateWindow()

    def change_faces_color_front(self, color):
        self.settings['faces.color:front'] = color
        self.view.update_vertex_buffer('faces.color:front', self.view.array_faces_color_front)
        self.view.updateGL()
        self.app.main.activateWindow()

    def change_faces_color_back(self, color):
        self.settings['faces.color:back'] = color
        self.view.update_vertex_buffer('faces.color:back', self.view.array_faces_color_back)
        self.view.updateGL()
        self.app.main.activateWindow()

    # ==========================================================================
    # camera
    # ==========================================================================

    # ==========================================================================
    # tools
    # ==========================================================================

    # ==========================================================================
    # my assembly model
    # ==========================================================================

    def create_floor_assembly(self, mu=0.4, cm=0, delta=0, massopt=False):
        self.assembly = create_assembly(mu, cm, delta, massopt)
        self.view.glInit()
        self.view.make_buffers()
        self.view.updateGL()
        self.app.main.activateWindow()


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

    pass
