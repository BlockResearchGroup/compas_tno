from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas_tessellation.app.core import App

from compas_tessellation.app.app1.view import View
from compas_tessellation.app.app1.controller import Controller

from compas_tessellation.app.app1 import CONFIG
from compas_tessellation.app.app1 import STYLE

__all__ = ['AssemblyViewer']


class AssemblyViewer(App):
    """"""

    def __init__(self):
        super(AssemblyViewer, self).__init__(CONFIG, STYLE)
        self.controller = Controller(self)
        self.view = View(self.controller)
        self.setup()
        self.init()
        self.view.glInit()
        self.view.setup_grid()
        # self.view.setup_axes()
        self.controller.create_floor_assembly()

    @property
    def assembly(self):
        return self.controller.assembly

    @assembly.setter
    def assembly(self, assembly):
        self.controller.assembly = assembly
        # self.controller.center_assembly()
        self.view.glInit()
        self.view.make_buffers()
        self.view.updateGL()


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

    viewer = AssemblyViewer()
    viewer.show()
