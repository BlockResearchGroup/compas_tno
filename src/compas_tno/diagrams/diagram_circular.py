
import math

def create_circular_radial_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0, diagonal=False, partial_diagonal=False):
    """ Helper to construct a circular radial FormDiagram with hoops not equally spaced in plan.

    Parameters
    ----------
    center : list
        Planar coordinates of the form-diagram [xc, yc].

    radius: float
        Radius of the form-diagram

    discretisation : list
        Number of meridians, and of spikes from the center on the dome form-diagram.

    r_oculus: float
        Value of the radius of the oculus, if no oculus is present should be set to zero.

    diagonal: float
        Activate diagonal in the quads.

    partial_diagonal: float
        Activate partial diagonal in the quads.

    Returns
    -------
    obj
        FormDiagram.

    """

    xc = center[0]
    yc = center[1]
    n_radial = discretisation[0]
    n_spikes = discretisation[1]
    theta = 2*math.pi/n_spikes
    r_div = (radius - r_oculus)/n_radial
    lines = []

    for nr in range(n_radial+1):
        for nc in range(n_spikes):
            if (r_oculus + nr * r_div) > 0.0:
                # Meridian Elements
                xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                xb = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1))
                ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                yb = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1))
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            if nr <= n_radial - 1:
                # Radial Elements
                xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                xb = xc + (r_oculus + (nr + 1) * r_div) * math.cos(theta * nc)
                ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                yb = yc + (r_oculus + (nr + 1) * r_div) * math.sin(theta * nc)
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

    if diagonal:

        for nr in range(n_radial):
            for nc in range(n_spikes):
                if (r_oculus + nr * r_div) > 0.0:

                    # Meridian Element i
                    xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                    xb = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1))
                    ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                    yb = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1))

                    # Meridian Element i + 1
                    xa_ = xc + (r_oculus + (nr + 1) * r_div) * math.cos(theta * nc)
                    xb_ = xc + (r_oculus + (nr + 1) * r_div) * math.cos(theta * (nc + 1))
                    ya_ = yc + (r_oculus + (nr + 1) * r_div) * math.sin(theta * nc)
                    yb_ = yc + (r_oculus + (nr + 1) * r_div) * math.sin(theta * (nc + 1))


                    if partial_diagonal == True:
                        if nc + 1 > n_spikes/2:
                            lines.append([[xa, ya, 0.0], [xb_, yb_, 0.0]])
                        else:
                            lines.append([[xa_, ya_, 0.0], [xb, yb, 0.0]])
                    else:
                        lines.append([[xa, ya, 0.0], [xb_, yb_, 0.0]])
                        lines.append([[xa_, ya_, 0.0], [xb, yb, 0.0]])

    form = cls.from_lines(lines, delete_boundary_edges=True)

    if r_oculus:
        for key in form.faces():
            centroid = form.face_centroid(key)
            if centroid[0] == xc and centroid[1] == yc:
                form.set_face_attribute(key, 'is_loaded', False)

    [bnds] = form.vertices_on_boundaries()
    for key in bnds:
        form.vertex_attribute(key, 'is_fixed', True)

    for u,v in form.edges_on_boundary():
        form.edge_attribute((u,v), 'is_edge', False)

    return form


def create_circular_radial_spaced_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0, diagonal=False, partial_diagonal=False):
    """ Helper to construct a circular radial FormDiagram with hoops not equally spaced in plan, but equally spaced with regards to the projection on a hemisphere.

    Parameters
    ----------
    center : list
        Planar coordinates of the form-diagram [xc, yc].

    radius: float
        Radius of the form-diagram

    discretisation : list
        Number of meridians, and of spikes from the center on the dome form-diagram.

    r_oculus: float
        Value of the radius of the oculus, if no oculus is present should be set to zero.

    diagonal: float
        Activate diagonal in the quads.

    partial_diagonal: float
        Activate partial diagonal in the quads.

    Returns
    -------
    obj
        FormDiagram.

    """

    xc = center[0]
    yc = center[1]
    n_radial = discretisation[0]
    n_spikes = discretisation[1]
    theta = 2*math.pi/n_spikes
    r_div = (radius - r_oculus)/n_radial
    lines = []

    for nr in range(n_radial+1):
        for nc in range(n_spikes):
            if (r_oculus + nr) > 0:
                # Meridian Elements
                xa = xc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.cos(theta * nc)
                xb = xc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.cos(theta * (nc + 1))
                ya = yc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.sin(theta * nc)
                yb = yc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.sin(theta * (nc + 1))
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            if nr <= n_radial - 1:

                # Radial Elements
                xa = xc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.cos(theta * nc)
                xb = xc + (r_oculus + radius * math.cos((n_radial - (nr + 1))/n_radial * math.pi/2)) * math.cos(theta * nc)
                ya = yc + (r_oculus + radius * math.cos((n_radial - nr)/n_radial * math.pi/2)) * math.sin(theta * nc)
                yb = yc + (r_oculus + radius * math.cos((n_radial - (nr + 1))/n_radial * math.pi/2)) * math.sin(theta * nc)
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

    if diagonal:

        for nr in range(n_radial):
            for nc in range(n_spikes):
                if (r_oculus + nr * r_div) > 0.0:

                    # Meridian Element i
                    xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                    xb = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1))
                    ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                    yb = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1))

                    # Meridian Element i + 1
                    xa_ = xc + (r_oculus + (nr + 1) * r_div) * math.cos(theta * nc)
                    xb_ = xc + (r_oculus + (nr + 1) * r_div) * math.cos(theta * (nc + 1))
                    ya_ = yc + (r_oculus + (nr + 1) * r_div) * math.sin(theta * nc)
                    yb_ = yc + (r_oculus + (nr + 1) * r_div) * math.sin(theta * (nc + 1))

                    if partial_diagonal == True:
                        if nc + 1 > n_spikes/2:
                            lines.append([[xa, ya, 0.0], [xb_, yb_, 0.0]])
                        else:
                            lines.append([[xa_, ya_, 0.0], [xb, yb, 0.0]])
                    else:
                        lines.append([[xa, ya, 0.0], [xb_, yb_, 0.0]])
                        lines.append([[xa_, ya_, 0.0], [xb, yb, 0.0]])

    form = cls.from_lines(lines, delete_boundary_edges=True)

    if r_oculus:
        for key in form.faces():
            centroid = form.face_centroid(key)
            if centroid[0] == xc and centroid[1] == yc:
                form.set_face_attribute(key, 'is_loaded', False)

    [bnds] = form.vertices_on_boundaries()
    for key in bnds:
        form.vertex_attribute(key, 'is_fixed', True)

    return form


def create_circular_spiral_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0):
    """ Helper to create a spiral Form-diagram and set common attributes.

    Parameters
    ----------
    center : list
        Planar coordinates of the form-diagram [xc, yc].

    radius: float
        Radius of the form-diagram

    discretisation : list
        Number of meridians, and of spikes from the center on the dome form-diagram.

    r_oculus: float
        Value of the radius of the oculus, if no oculus is present should be set to zero.

    Returns
    -------
    obj
        FormDiagram.

    """

    xc = center[0]
    yc = center[1]
    n_radial = discretisation[0]
    n_spikes = discretisation[1]
    theta = 2*math.pi/n_spikes
    r_div = (radius - r_oculus)/n_radial
    lines = []

    for nr in range(n_radial + 1):
        for nc in range(n_spikes):
            if nr > 0.0:  # This avoid the center...
                if nr % 2 == 0:
                    # Diagonal to Up
                    xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                    xb = xc + (r_oculus + (nr - 1) * r_div) * math.cos(theta * (nc + 1/2))
                    ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                    yb = yc + (r_oculus + (nr - 1) * r_div) * math.sin(theta * (nc + 1/2))

                    # Diagonal to Down
                    xa_ = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                    xb_ = xc + (r_oculus + (nr - 1) * r_div) * math.cos(theta * (nc - 1/2))
                    ya_ = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                    yb_ = yc + (r_oculus + (nr - 1) * r_div) * math.sin(theta * (nc - 1/2))

                    lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
                    lines.append([[xa_, ya_, 0.0], [xb_, yb_, 0.0]])
                else:
                    # Diagonal to Up
                    xa = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1/2))
                    xb = xc + (r_oculus + (nr - 1) * r_div) * math.cos(theta * (nc + 1))
                    ya = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1/2))
                    yb = yc + (r_oculus + (nr - 1) * r_div) * math.sin(theta * (nc + 1))

                    # Diagonal to Down
                    xa_ = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1/2))
                    xb_ = xc + (r_oculus + (nr - 1) * r_div) * math.cos(theta * (nc))
                    ya_ = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1/2))
                    yb_ = yc + (r_oculus + (nr - 1) * r_div) * math.sin(theta * (nc))

                    lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
                    lines.append([[xa_, ya_, 0.0], [xb_, yb_, 0.0]])
                if nr == n_radial:
                    xa = xc + (r_oculus + nr * r_div) * math.cos(theta * nc)
                    xb = xc + (r_oculus + nr * r_div) * math.cos(theta * (nc + 1))
                    ya = yc + (r_oculus + nr * r_div) * math.sin(theta * nc)
                    yb = yc + (r_oculus + nr * r_div) * math.sin(theta * (nc + 1))
                    lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
            if nr == 0 and r_oculus > 0.0:
                # If oculus, this will be the compression ring
                xa = xc + (r_oculus) * math.cos(theta * (nc))
                xb = xc + (r_oculus) * math.cos(theta * (nc + 1))
                ya = yc + (r_oculus) * math.sin(theta * (nc))
                yb = yc + (r_oculus) * math.sin(theta * (nc + 1))
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

    form = cls.from_lines(lines, delete_boundary_edges=True)

    if r_oculus:
        for key in form.faces():
            centroid = form.face_centroid(key)
            if centroid[0] == xc and centroid[1] == yc:
                form.set_face_attribute(key, 'is_loaded', False)

    [bnds] = form.vertices_on_boundaries()
    for key in bnds:
        form.vertex_attribute(key, 'is_fixed', True)

    form.delete_boundary_edges() # Check what happens if there is oculus


    return form