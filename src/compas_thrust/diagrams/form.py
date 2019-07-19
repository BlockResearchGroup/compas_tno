
from compas.utilities import geometric_key
from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import vertical_from_zmax
from compas_thrust.plotters.plotters import plot_form
from compas_thrust.plotters.plotters import plot_force
from random import shuffle

def _form(form, keep_q=False):

    """ s the FormDiagram by shuffling the edges.

    Parameters
    ----------
    form : obj
        Original FormDiagram.

    Returns
    -------
    obj
        Shuffled FormDiagram.

    """

    # Edges

    edges = [form.edge_coordinates(u, v) for u, v in form.edges()]
    edges = [[sp[:2] + [0], ep[:2] + [0]] for sp, ep in edges]
    qs = {geometric_key(form.edge_midpoint(u, v)[:2] + [0]) : form.get_edge_attribute((u,v), 'q') for u, v in form.edges()}
    shuffle(edges)

    form_ = FormDiagram.from_lines(edges, delete_boundary_face=False)
    form_.update_default_edge_attributes({'is_symmetry': False})
    sym = [geometric_key(form.edge_midpoint(u, v)[:2] + [0])for u, v in form.edges_where({'is_symmetry': True})]
    for u, v in form_.edges():
        if geometric_key(form_.edge_midpoint(u, v)) in sym:
            form_.set_edge_attribute((u, v), 'is_symmetry', True)
        if keep_q:
            form_.set_edge_attribute((u, v), 'q', qs[geometric_key(form.edge_midpoint(u, v)[:2] + [0])])



    # Vertices

    gkey_key = form_.gkey_key()
    for key, vertex in form.vertex.items():
        gkey = geometric_key(form.vertex_coordinates(key)[:2] + [0])
        form_.vertex[gkey_key[gkey]] = vertex

    form_.attributes['indset'] = []

    return form_

def adapt_tna(form, zmax = 5.0, plot = False, delete_face = False):

    if delete_face:
        form.delete_face(0)
    corners = list(form.vertices_where({'is_fixed': True}))
    form.set_vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
    form.update_boundaries(feet=2)
    form.plot()
    force  = ForceDiagram.from_formdiagram(form)
    horizontal(form, force, alpha=100, display=False)
    vertical_from_zmax(form,zmax,display=False)

    if plot:
        plot_form(form).show()
        plot_force(force).show()

    return form

def evaluate_a(form, plot=True):

    a_total = 0
    a_max = 0
    for u, v, attr in form.edges_where({'is_edge': True}, True):
        a = attr['a']
        a_total += a
        l = form.edge_length(u,v)
        a = a*l
        if a > a_max:
            a_max = a
    if plot is True:
        print('Angle Deviation  Max: {0}'.format(a_max))

    return a_max
