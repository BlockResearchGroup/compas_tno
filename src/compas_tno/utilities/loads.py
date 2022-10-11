
from compas_tno.shapes.dome import dome_zt_update
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update
from compas_tno.utilities.interpolation import get_shape_middle_pattern


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


def apply_selfweight_from_shape(form, shape, pz_negative=True, normalize=True):
    """Apply selfweight to the nodes of the form diagram based on the shape

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the selfweight
    shape : Shape
        Shape of the masonry
    pz_negative : bool, optional
        Wether or not the vertical loads are negative, by default True
    normalize : bool, optional
        Wether or not normalise the selfweight to match shape.total_weight, by default True

    Returns
    -------
    None
        The FormDiagram is modified in place
    """

    form_ = form.copy()
    total_selfweight = shape.compute_selfweight()
    ro = shape.ro
    thk = shape.datashape['thk']

    x = form.vertices_attribute('x')
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

    factor = 1.0 * ro * thk  # Transform tributary area in tributary load
    if normalize:
        factor = total_selfweight/pzt
    if pz_negative:
        factor *= -1  # make loads negative

    for key in form.vertices():
        pzi = factor * form.vertex_attribute(key, 'pz')
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
        # if pz_negative:
        #     pz *= -1  # make loads negative
        form.vertex_attribute(key, 'pz', value=pz)
        pzt += pz

    if plot:
        print('total load applied:', pzt)


def apply_selfweight_from_thrust(form, thickness=0.5, density=20.0):
    """Lump the selfweight in the nodes of the thrust network based on their current position.
    The loads are computed based on the tributary area times the thickness times the density.
    For variable thickness the nodal attribute `thk` is considered.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to be considered
    thickness : float, optional
        The thickness of the problem, by default 0.50
        If None is passed, the thickness is taken from the nodal attribute `thk`
    density : float, optional
        The density of the material, by default 20.0

    Return
    ------
    None
        The form diagram is updated in place
    """

    for key in form.vertices():
        ai = form.vertex_area(key)
        if thickness:
            load = -1 * ai * thickness * density
        else:
            thk = form.vertex_attribute(key, 'thk')
            load = -1 * ai * thk * density
        form.vertex_attribute(key, 'pz', load)


def apply_horizontal_multiplier(form, lambd=1.0, direction='x'):
    """Modify the applied loads considering a load multiplier.

    Parameters
    ----------
    form : FormDiagram
        Form diagram to apply the horizontal multiplier.
    lambd : float, optional
        Value of the horizontal multiplier, by default ``1.0``.
    direction : str, optional
        Direction to apply the loads, ``x`` or ``y``, by default ``x``.

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

    raise NotImplementedError('Not implemented')
