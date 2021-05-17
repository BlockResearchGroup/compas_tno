from compas_tna.diagrams import ForceDiagram

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'ForceDiagram',
]


class ForceDiagram(ForceDiagram):

    """The ``Forcediagram`` class imports the attributes and set ups from compas_tna.diagrams.FormDiagram and include some functionalities useful for the assessment of masonry structures.

    It defined the form-diagram that will be the layout of the forces within the structure

    Notes
    -----
    A ``ForceDiagram`` is generated as the dual of the ``FormDiagram`` and it makes sense for structural analysis when both diagrams are reciprocal, i.e. when each edge in the Force diagram parallel to its corresponding edge in the Form diagram.

    """

    __module__ = 'compas_tna.diagrams'

    def __init__(self):
        super(ForceDiagram, self).__init__()

    # --------------------------------------------------------------------------
    # Convenience functions for retrieving attributes of the force diagram.
    # --------------------------------------------------------------------------

    def xy(self):
        """The XY coordinates of the vertices.

        Returns
        -------
        list
        """
        return self.vertices_attributes('xy')

    def fixed(self):
        """The identifiers of the fixed vertices.

        Returns
        -------
        list
        """
        return list(self.vertices_where({'is_fixed': True}))


    def fixed_x(self):
        """The identifiers of the vertices fixed in ``x`` only.

        Returns
        -------
        list
        """
        return list(self.vertices_where({'is_fixed_x': True, 'is_fixed': False}))

    def fixed_y(self):
        """The identifiers of the vertices fixed in ``y`` only.

        Returns
        -------
        list
        """
        return list(self.vertices_where({'is_fixed_y': True, 'is_fixed': False}))


    def anchor(self):
        """Get an anchor to the force diagram.

        Returns
        -------
        int
        """
        return next(self.vertices())


# def update_forcediagram(form,force):

#     n = form.number_of_vertices()
#     print('Vertices on form {0}'.format(n))
#     k_i = form.key_index()
#     print(len(k_i))
#     xyz = zeros((n, 3))
#     for key in form.vertices():
#         i = k_i[key]
#         xyz[i, :] = form.vertex_coordinates(key)
#     xy = xyz[:, :2]

#     edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True})]
#     # edges = [[k_i[u], k_i[v]] for u, v in form.edges()]
#     C	 = connectivity_matrix(edges, 'csr')
#     edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True})]
#     # q = [attr['q'] for u, v, attr in form.edges(True)]
#     q = [form.edge_attribute((u,v),'q') for u, v in form.edges_where({'_is_edge': True})]
#     Q = diags(q)
#     uv = C.dot(xy)

#     _k_i = force.key_index()
#     _known = []
#     for x in force.fixed():
#         _known.append(_k_i[x])

#     _n = force.number_of_vertices()
#     _xyz = zeros((_n, 3))
#     for key in force.vertices():
#         i = _k_i[key]
#         _xyz[i, :] = force.vertex_coordinates(key)
#     _xy = _xyz[:, :2]

#     _edges = force.ordered_edges(form)
#     _C = connectivity_matrix(_edges, 'csr')
#     _Ct = _C.transpose()

#     _xy = spsolve_with_known(_Ct.dot(_C), _Ct.dot(Q).dot(uv), _xy, _known)

#     for key, attr in force.vertices(True):
#         i = _k_i[key]
#         attr['x'] = _xy[i, 0]
#         attr['y'] = _xy[i, 1]

#     return force

# def recalculate_qs(form,force):

#     k_i	 = form.key_index()
#     uv_i	= form.uv_index()
#     vcount  = form.number_of_vertices()
#     anchors = list(form.anchors())
#     fixed   = list(form.fixed())
#     fixed   = set(anchors + fixed)
#     fixed   = [k_i[key] for key in fixed]
#     free	= list(set(range(vcount)) - set(fixed))
#     edges   = [(k_i[u], k_i[v]) for u, v in form.edges()]
#     xyz	 = array(form.get_vertices_attributes('xyz'), dtype='float64')
#     for i in range(vcount):
#         xyz[i,2] = 0.0
#     C	   = connectivity_matrix(edges, 'csr')
#     _scale = force.scale
#     _xyz   = array(force.get_vertices_attributes('xyz'), dtype='float64')
#     for i in range(_xyz.shape[1]):
#         _xyz[i,2] = 0.0
#     _edges = force.ordered_edges(form)
#     _C	 = connectivity_matrix(_edges, 'csr')
#     uvw  = C.dot(xyz)
#     _uvw = _C.dot(_xyz)
#     l	= normrow(uvw)
#     _l   = normrow(_uvw)
#     f	= _scale * _l
#     q	= _l / l

#     print('q from force')

#     return q
