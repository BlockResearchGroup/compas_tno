
from compas_tno.shapes.dome import dome_zt_update
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update
from compas_tno.utilities.interpolation import get_shape_middle_pattern


__all__ = [
    'apply_selfweight_from_shape',
    'apply_selfweight_from_pattern',
    'apply_horizontal_multiplier',
    'apply_fill_load',

    'apply_selfweight_from_shape_proxy',
]


def apply_selfweight_from_shape_proxy(formdata, shapedata):
    """Apply selfweight to the nodes of the form diagram based on the shape from proxy

    Parameters
    ----------
    formdata : dict
        Data of the form diagram to apply the selfweight.
    shapedata : Shape
        Data of the shape of the masonry.

    Returns
    -------
    formdata
        Data of the form diagram with loads applied
    """

    from compas_tno.diagrams import FormDiagram
    from compas_tno.shapes import Shape

    form = FormDiagram.from_data(formdata)
    shape = Shape.from_data(shapedata)

    apply_selfweight_from_shape(form, shape)

    return form.to_data()


def apply_selfweight_from_shape(form, shape, pz_negative=True):
    """Apply selfweight to the nodes of the form diagram based on the shape

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the selfweight
    shape : Shape
        Shape of the masonry
    pz_negative : bool, optional
        Wether or not the vertical loads are negative, by default True

    Returns
    -------
    None
        The FormDiagram is modified in place
    """

    form_ = form.copy()
    total_selfweight = shape.compute_selfweight()

    x = form.vertices_attribute('x')  # check if array is necessary here
    y = form.vertices_attribute('y')

    if shape.datashape['type'] == 'dome':
        zt = dome_zt_update(x, y, shape.datashape['radius'], shape.datashape['t'], shape.datashape['center'])
    elif shape.datashape['type'] == 'crossvault':
        zt = crossvault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        zt = pointed_vault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'],
                                         hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
    else:
        XY = form.vertices_attributes('xy')
        zt = get_shape_middle_pattern(shape, XY)

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


def apply_selfweight_from_pattern(form, pattern, plot=False, pz_negative=True, tol=10e-4):
    """Apply selfweight to the nodes considering a different Form Diagram to locate loads.
    Note: the base pattern has to coincide with nodes from the original form diagram.

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the selfweight.
    pattern : Mesh
        Mesh in which the forces should be based.
    plot : bool, optional
        Whether or not plot the diagrams and its overlap, by False
    pz_negative : bool, optional
        Wether or not the vertical loads are negative by default ``True``.
    tol : float, optional
        Tolerance for the nodal position match of form and pattern, by default ``10e-4``.

    Returns
    -------
    None
        The FormDiagram is modified in place.
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

    if plot:
        print('total load applied:', pzt)


def apply_horizontal_multiplier(form, lambd=0.1, direction='x'):
    """Modify the applied loads considering a load multiplier.

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the horizontal multiplier.
    lambd : float, optional
        Value of the horizontal multiplier, by default ``0.1``.
    direction : str, optional
        direction to apply the loads, by default ``x``.

    Returns
    -------
    None
        The FormDiagram is modified in place.
    """

    arg = 'p' + direction

    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        form.vertex_attribute(key, arg, -1 * pz * lambd)  # considers that swt (pz) is negative


def apply_fill_load(form):
    """Modify the applied loads considering a fill.

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the horizontal multiplier.

    Note
    -------
        In development.
    """

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
