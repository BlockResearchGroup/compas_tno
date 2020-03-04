"""Cablenet equilibrium.

- make a network from sample data
- set default vertex and edge attributes
- identify *anchored* vertices
- convert network data to numerical data
- use a *key-index map* to reference items in the numerical data
- run the force density method
- update the network
- make a viewer and display the result


Note
----
This examples requires PyOpenGL for visualization.


"""

import compas

from compas.datastructures import Network
from compas_viewers.meshviewer import MeshViewer
from compas_tno.plotters import plot_form
from compas_tna.diagrams import FormDiagram

from compas.numerical import fd_numpy

file = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_TNA.json'
form = FormDiagram.from_json(file)
plot_form(form)

viewer = MeshViewer()
viewer.mesh = form
# viewer.setup()
viewer.show()
