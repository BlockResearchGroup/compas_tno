from compas_tna.diagrams import ForceDiagram
from compas.datastructures import mesh_dual


class ForceDiagram(ForceDiagram):
    """The ``ForceDiagram`` represents the equilibrium of the forces in the form diagram.

    Parameters
    ----------


    Attributes
    ----------

    """

    def __init__(self):
        super(ForceDiagram, self).__init__()
        self.attributes.update({
            'name': 'Force'})
        self.update_default_vertex_attributes({
            'is_fixed': False,
            'is_param': False})
        self.update_default_edge_attributes({
            'l': 0.0,
            'target_vector': None,
            'target_length': None})

    # --------------------------------------------------------------------------
    # Constructors
    # --------------------------------------------------------------------------

    @classmethod
    def from_formdiagram(cls, formdiagram):
        """Construct a force diagram from a form diagram.

        Parameters
        ----------
        formdiagram : :class:`~compas_tno.diagrams.ForceDiagram`
            The form diagram.

        Returns
        -------
        :class:`~compas_tno.diagrams.ForceDiagram`
        """
        forcediagram = mesh_dual(formdiagram, cls)
        forcediagram.dual = formdiagram
        formdiagram.dual = forcediagram
        return forcediagram

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

    def anchor(self):
        """Get an anchor to the force diagram.

        Returns
        -------
        int
        """
        return next(self.vertices())

    # --------------------------------------------------------------------------
    # Plot
    # --------------------------------------------------------------------------

    def plot(self):
        """Plot a force diagram with a plotter with all the default settings."""
        from compas_plotters import Plotter
        plotter = Plotter(figsize=(8, 8))
        plotter.add(self)
        plotter.show()

    # --------------------------------------------------------------------------
    # Helpers
    # --------------------------------------------------------------------------

    def edges_where_dual(self, conditions, data=False):
        """Get edges for which a certain condition or set of conditions is true for the corresponding edges in the diagram's dual.

        Parameters
        ----------
        conditions : dict
            A set of conditions in the form of key-value pairs.
            The keys should be attribute names. The values can be attribute
            values or ranges of attribute values in the form of min/max pairs.
        data : bool, optional
            Yield the edges and their data attributes.
            Default is ``False``.

        Yields
        ------
        2-tuple
            The next edge as a (u, v) tuple, if ``data=False``.
            The next edge as a ((u, v), data) tuple, if ``data=True``.

        Examples
        --------
        >>>

        """
        for edge in list(self.edges()):
            is_match = True

            dual_edge = self.dual_edge(edge)
            dual_edge_attr = self.dual.edge_attributes(dual_edge)

            for cond_name, cond_value in conditions.items():
                dual_method = getattr(self.dual, cond_name, None)

                if dual_method and callable(dual_method):
                    dual_value = dual_method(dual_edge)
                elif cond_name in dual_edge_attr:
                    dual_value = dual_edge_attr[cond_name]
                else:
                    is_match = False
                    break

                if isinstance(dual_value, list):
                    if cond_value not in dual_value:
                        is_match = False
                        break
                elif isinstance(cond_value, (tuple, list)):
                    minval, maxval = cond_value
                    if dual_value < minval or dual_value > maxval:
                        is_match = False
                        break
                else:
                    if cond_value != dual_value:
                        is_match = False
                        break

            if is_match:
                if data:
                    yield edge, self.edge_attributes(edge)
                else:
                    yield edge

    def dual_edge(self, edge):
        """Find the cooresponding edge in the diagram's dual.

        Parameters
        ----------
        edge : tuple of int
            The edge identifier.

        Returns
        -------
        tuple (int, int) or None
            The identifier of the dual edge if it exists.
        """
        for u, v in self.dual.face_halfedges(edge[0]):
            if self.dual.halfedge[v][u] == edge[1]:
                if self.dual.has_edge((u, v)):
                    return u, v
                return v, u

    def is_dual_edge_external(self, edge):
        """Verify if the corresponding edge in the diagram's dual is marked as "external".

        Parameters
        ----------
        edge : tuple of int
            The edge identifier.

        Returns
        -------
        bool
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'is_external')

    def is_dual_edge_reaction(self, edge):
        """Verify if the corresponding edge in the diagram's dual is marked as "reaction".

        Parameters
        ----------
        edge : tuple of int
            The edge identifier.

        Returns
        -------
        bool
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'is_reaction')

    def is_dual_edge_load(self, edge):
        """Verify if the corresponding edge in the diagram's dual is marked as "load".

        Parameters
        ----------
        edge : tuple of int
            The edge identifier.

        Returns
        -------
        bool
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'is_load')

    def is_dual_edge_ind(self, edge):
        """Verify if the corresponding edge in the diagram's dual is marked as "independent".

        Parameters
        ----------
        edge : tuple of int
            The edge identifier.

        Returns
        -------
        bool
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'is_ind')

    def dual_edge_force(self, edge):
        """Retrieve the force in the corresponding edge of the diagram's dual.

        Parameters
        ----------
        edge : tuple(int, int)
            The edge identifier.

        Returns
        -------
        float
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'f')

    def dual_edge_angledeviation(self, edge):
        """Retrieve the angle deviation in the corresponding edge of the diagram's dual.

        Parameters
        ----------
        edge : tuple(int, int)
            The edge identifier.

        Returns
        -------
        float
        """
        return self.dual.edge_attribute(self.dual_edge(edge), 'a')

    def edge_index(self, form=None):
        """Construct a mapping between the identifiers of edges and the corresponding indices in a list of edges.

        Parameters
        ----------
        form : :class:`compas_ags.diagrams.FormDiagram`, optional
            If a form diagra is provided as reference, the list of edges is ordered such that it corresponds
            to the natural ordering of edges in the form diagram.

        Returns
        -------
        dict
            Mapping between edge identifiers and the correpsonding indices of the edges in a list.

        """
        if not form:
            return {edge: index for index, edge in enumerate(self.edges())}
            # return {edge: index for index, edge in enumerate(form.edges_where({'_is_edge': True}))}
        edge_index = dict()
        # for index, (u, v) in enumerate(form.edges()):
        for index, (u, v) in enumerate(form.edges_where({'_is_edge': True})):
            f1 = form.halfedge[u][v]
            f2 = form.halfedge[v][u]
            edge_index[f1, f2] = index
            # the weird side-effect of this is that edges get rotated if necessary
        return edge_index

    def ordered_edges(self, form):
        """"Construct a list of edges with the same order as the corresponding edges of the form diagram.

        Parameters
        ----------
        form : :class:`compas_ags.diagrams.FormDiagram`

        Returns
        -------
        list
        """
        edge_index = self.edge_index(form=form)
        index_edge = {index: edge for edge, index in edge_index.items()}
        edges = [index_edge[index] for index in range(self.number_of_edges())]
        return edges
