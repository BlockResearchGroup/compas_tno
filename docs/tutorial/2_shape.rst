.. _shape:

********************************************************************************
Shape
********************************************************************************

.. currentmodule:: compas_tno.shapes

.. currentmodule:: compas_tno.shapes.Shape

.. highlight:: python

This tutorial provides a quick tour of the generation of :mod:`Shape <compas_tno.shapes.Shape>`.

The Shape can be created through a series of methods, such as:

* :mod:`from_meshes <compas_tno.shapes.Shape.from_meshes>`: create the Shape by providing meshes representing the intrados and extrados.
* :mod:`from_pointcloud <compas_tno.shapes.Shape.from_pointcloud>`: create the Shape from pointclouds.
* :mod:`from_library <compas_tno.shapes.Shape.from_library>`: which creates a series of parametric Shapes that correspond to common masonry vaulted geometries.

The TNO library of Shapes is based on commom layouts which will be described herein.

Shapes Library
--------------------

.. figure:: ../_images/shapes.png
    :figclass: figure
    :class: figure-img img-fluid

The methods to create each of these shapes is described herein.

a) :mod:`create_dome <compas_tno.shapes.Shape.create_dome>`: hemispheric dome geometry defined from its radius :math:`R`, thickness :math:`t`, center position :math:`[x_\mathrm{c}, y_\mathrm{c}]`, etc. One example is presented below:

.. code-block:: Python

    from compas_tno.shapes import Shape
    shape = Shape.create_dome(center=[5.0, 5.0],
                              radius=5.0,
                              thk=0.5,
                              r_oculus=0.75)

b) :mod:`create_crossvault <compas_tno.shapes.Shape.create_crossvault>`: rounded cross vault geometry defined from the span :math:`s`, thickness :math:`t`, springing angle :math:`\beta`, corner positions :math:`[[x_\mathrm{0}, y_\mathrm{0}], [x_\mathrm{f}, y_\mathrm{f}]]`. One example is presented below:

.. code-block:: Python

    from compas_tno.shapes import Shape
    shape = Shape.create_crossvault(xy_span=[[0.0, 10.0], [0.0, 10.0]],
                                    thk=0.5,
                                    spr_angle=30.0)

c) :mod:`create_pointedcrossvault <compas_tno.shapes.Shape.create_pointedcrossvault>`: pointed cross vault geometry defined from the span :math:`s`, thickness :math:`t`, springing angle :math:`\beta`, corner positions :math:`[[x_\mathrm{0}, y_\mathrm{0}], [x_\mathrm{f}, y_\mathrm{f}]]`, height of center point :math:`h_\mathrm{c}`, and height of the boundary points :math:`h_\mathrm{b}`. One example is presented below:

.. code-block:: Python

    from compas_tno.shapes import Shape
    shape = Shape.create_pointedcrossvault(xy_span=[[0.0, 10.0], [0.0, 10.0]],
                                           thk=0.5,
                                           hc=7.0)

d) :mod:`create_pavillionvault <compas_tno.shapes.Shape.create_pavillionvault>`: rounded pavillion vault geometry defined from the span :math:`s`, thickness :math:`t`, springing angle :math:`\beta`, corner positions :math:`[[x_\mathrm{0}, y_\mathrm{0}], [x_\mathrm{f}, y_\mathrm{f}]]`. One example is presented below:

.. code-block:: Python

    from compas_tno.shapes import Shape
    shape = Shape.create_pavillionvault(xy_span=[[0.0, 10.0], [0.0, 10.0]],
                                        thk=0.5,
                                        spr_angle=30.0)

To further explore the library of parametric shapes, check the :mod:`Shapes <compas_tno.shapes.Shape>` full documentation.

TNO Viewer
-------------

For 3D visualisation, the :mod:`Viewer <compas_tno.viewers.TNOViewer>` can be used. This `COMPAS View2 <https://compas.dev/compas_view2/latest/>`_ based visualisation enable to add multiple 3D objects in the scene, such as the shape, the thrust network, cracks, reaction forces, loads, etc. The viewer is used in the :ref:`examples <examples>` of this tutorial. The code to visualise the shape of a crossvault is:

.. code-block:: Python

    from compas_tno.viewers import Viewer
    from compas_tno.shapes import Shape
    vault = Shape.create_crossvault(xy_span=[[0.0, 10.0], [0.0, 10.0]],
                                    thk=0.5,
                                    spr_angle=30.0)
    view = View(shape=vault)
    view.draw_shape()
    view.show()
