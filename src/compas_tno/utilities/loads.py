
from compas_tno.shapes.dome import dome_zt_update
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update


__all__ = [
    'apply_selfweight_from_shape',
    'apply_selfweight_from_pattern',
    'apply_horizontal_multiplier',
    'apply_fill_load',

    'apply_selfweight_from_shape_proxy',
]


def apply_selfweight_from_shape_proxy(formdata, shapedata):  # this works only forlibrary shapes
    # TODO: crate a proper to_data from_data for shapes, and make it happen.

    from compas_tno.diagrams import FormDiagram
    from compas_tno.shapes import Shape

    form = FormDiagram.from_data(formdata)
    shape = Shape.from_library(shapedata)

    apply_selfweight_from_shape(form, shape)

    return form.to_data()


def apply_selfweight_from_shape(form, shape, pz_negative=True):
    """Apply selfweight to the nodes of the form diagram based on the shape"""

    form_ = form.copy()
    total_selfweight = shape.compute_selfweight()

    x = form.vertices_attribute('x')  # check if array is necessary here
    y = form.vertices_attribute('y')

    if shape.datashape['type'] == 'dome':
        zt = dome_zt_update(x, y, shape.datashape['radius'], shape.datashape['t'], shape.datashape['center'])
    elif shape.datashape['type'] == 'crossvault':
        zt = crossvault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        zt = pointed_vault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'], hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
    else:
        XY = form.vertices_attributes('xy')
        zt = shape.get_middle_pattern(XY)

    i = 0
    for key in form_.vertices():
        z = float(zt[i])
        form_.vertex_attribute(key, 'z', value=z)
        form.vertex_attribute(key, 'target', value=z)
        i += 1

    pzt = 0
    for key in form.vertices():
        pz = form_.vertex_area(key)
        form.vertex_attribute(key, 'pz', value=pz)
        pzt += pz

    if shape.datashape['type'] == 'arch' or shape.datashape['type'] == 'pointed_arch':
        pzt = 0
        for key in form.vertices():
            form.vertex_attribute(key, 'pz', value=1.0)
            if form.vertex_attribute(key, 'is_fixed') is True:
                form.vertex_attribute(key, 'pz', value=0.5)
            pzt += form.vertex_attribute(key, 'pz')

    factor = total_selfweight/pzt

    for key in form.vertices():
        pzi = factor * form.vertex_attribute(key, 'pz')
        if pz_negative:
            pzi *= -1  # make loads negative
        form.vertex_attribute(key, 'pz', value=pzi)

    return


def apply_selfweight_from_pattern(form, pattern, plot=False, pz_negative=True, tol=10e-4):
    """Apply selfweight to the nodes considering a different Form Diagram to locate loads. Warning, the base pattern has to coincide with nodes from the original form diagram
    """

    form_ = pattern

    form.vertices_attribute('pz', 0.0)
    key_real_to_key = {}

    for key in form_.vertices():
        x, y, _ = form_.vertex_coordinates(key)
        for key_real in form.vertices():
            x_real, y_real, _ = form.vertex_coordinates(key_real)
            if x - tol < x_real < x + tol and y - tol < y_real < y + tol:
                key_real_to_key[key_real] = key
                break

    pzt = 0
    for key in key_real_to_key:
        pz = form_.vertex_attribute(key_real_to_key[key], 'pz')
        if pz_negative:
            pz *= -1  # make loads negative
        form.vertex_attribute(key, 'pz', value=pz)
        pzt += pz
    print('total load applied:', pzt)

    if plot:

        from compas_plotters import MeshPlotter

        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_vertices(text=key_real_to_key)
        plotter.show()

        plotter = MeshPlotter(form_, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_vertices(text={key: key for key in form_.vertices()})
        plotter.show()

        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 1) for key in form.vertices()})
        plotter.show()

        return


def apply_horizontal_multiplier(form, lambd=0.1, direction='x'):

    arg = 'p' + direction

    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        form.vertex_attribute(key, arg, -1 * pz * lambd)  # considers that swt (pz) is negative

    return


def apply_fill_load(form):

    print('Non implemented')

    return


# def vertex_projected_area(form, key):  # Modify to compute the projected aerea of all and save as an attribute
#     """Compute the projected tributary area of a vertex.

#     Parameters
#     ----------
#     key : int
#         The identifier of the vertex.

#     Returns
#     -------
#     float
#         The projected tributary area.

#     Example
#     -------
#     >>>

#     """

#     from compas.geometry import subtract_vectors
#     from compas.geometry import length_vector
#     from compas.geometry import cross_vectors

#     area = 0.

#     p0 = form.vertex_coordinates(key)
#     p0[2] = 0

#     for nbr in form.halfedge[key]:
#         p1 = form.vertex_coordinates(nbr)
#         p1[2] = 0
#         v1 = subtract_vectors(p1, p0)

#         fkey = form.halfedge[key][nbr]
#         if fkey is not None:
#             p2 = form.face_centroid(fkey)
#             p2[2] = 0
#             v2 = subtract_vectors(p2, p0)
#             area += length_vector(cross_vectors(v1, v2))

#         fkey = form.halfedge[nbr][key]
#         if fkey is not None:
#             p3 = form.face_centroid(fkey)
#             p3[2] = 0
#             v3 = subtract_vectors(p3, p0)
#             area += length_vector(cross_vectors(v1, v3))

#     return 0.25 * area
