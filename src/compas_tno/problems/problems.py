import time

# from numpy import asarray
from dataclasses import dataclass
from dataclasses import field
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

import numpy.typing as npt
from numpy import array
from numpy import hstack
from numpy import identity
from numpy import vstack
from numpy import zeros
from numpy.linalg import pinv

# from numpy.linalg import svd
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import vstack as svstack

from compas.geometry import Point
from compas.matrices import connectivity_matrix
from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import Envelope

# from compas_tno.algorithms import check_independents
from compas_tno.algorithms import check_horizontal_loads
from compas_tno.algorithms import find_independents
from compas_tno.utilities import apply_radial_symmetry
from compas_tno.utilities import apply_symmetry_from_axis
from compas_tno.utilities import build_symmetry_transformation
from compas_tno.utilities import find_sym_axis_in_rect_patterns


@dataclass
class Problem:
    """
    The ``Problem`` class stores the matrices used in the optimisation. These are listed as parameters of the class and described below.

    Parameters
    ----------

    None

    Attributes
    ----------

    q : array(m x 1)
        The vector of force densities
    m : int
        The number of edges
    n : int
        The number of vertices
    nb : int
        The number of fixed vertices
    E : array(2n x m)
        The horizontal equilibrium matrix
    C : array(n x m)
        The connectivity matrix
    Ct : array(m x n)
        The transposed connectivity matrix
    Ci : array(m x ni)
        The sliced connectivity matrix with regards to the internal vertices
    Cit : array(ni x m)
        The transpose of Ci
    Cb : array(nb x m)
        The sliced connectivity matrix with regards to the constrained (support) vertices
    U : array(m x m)
        The diagonal matrix of coorrdinate differences in the x-direction
    V : array(m x m)
        The diagonal matrix of coorrdinate differences in the y-direction
    P : array(n x 3)
        The applied external loads in all nodes
    free : list
        The list with the index of the free vertices
    fixed : list
        The list with the index of the fixed vertices
    phfree : array(ni x 2)
        The applied horizonal loads on the free vertices
    ph : array(ni x 2)
        A column vector with the the applied horizonal loads on the free vertices
    lb : array(n x 1)
        The lower-bound (intrados) limit for the nodes
    ub : array(n x 1)
        The upper-bound (extrados) limit for the nodes
    lb0 : array(n x 1)
        The lower-bound (intrados) limit for the nodes at the starting point (can be modified if thickness in the objective)
    ub0 : array(n x 1)
        The upper-bound (extrados) limit for the nodes at the starting point (can be modified if thickness in the objective)
    s : array(n x 1)
        The middle surface of the nodes
    X : array(n x 3)
        The nodal position of the vertices of the network
    x0 : array(n x 1)
        The x-position of the vertices in the network
    y0 : array(n x 1)
        The y-position of the vertices in the network
    free_x : list
        index of the vertices free to move in x
    free_y : list
        index of the vertices free to move in y
    rol_x : list
        index of the vertices constrained partially on x.
    rol_y : list
        index of the vertices constrained partially on y.
    Citx : array
        Connectivity matrix sliced transposed on the vertices free on x.
    City : array
        Connectivity matrix sliced transposed on the vertices free on y.
    Cbtx : array
        Connectivity matrix sliced transposed on the vertices partially fixed on x.
    Cbty : array
        Connectivity matrix sliced transposed on the vertices partially fixed on y.
    xlimits : array(n x 1)
        Limits on the x-direction in which the nodes can move
    ylimits : array(n x 1)
        Limits on the y-direction in which the nodes can move
    qmin : array(m x 1)
        Lower-bounds of the force densities in the edges
    qmax : array(m x 1)
        Upper-bounds of the force densities in the edges
    k_i : dict
        key-index dictionary
    uv_i : dict
        uv-index dictionary
    i_uv : dict
        index-uv dictionary
    ind : list
        List with the index of the independent edges
    k : int
        Number of independents in the problem
    dep : list
        List with the index of the dependent edges
    B : array(m x k)
        Matrix transforming the force densities in the independent edges to all force densities
    rho : float
        Density of the material
    constraints : list
        List with the constraints of the problem

    """

    q: Optional[npt.NDArray] = None
    m: Optional[int] = None
    n: Optional[int] = None
    ni: Optional[int] = None
    nb: Optional[int] = None
    E: Optional[npt.NDArray] = None
    C: Optional[npt.NDArray] = None
    Ct: Optional[npt.NDArray] = None
    Ci: Optional[npt.NDArray] = None
    Cit: Optional[npt.NDArray] = None
    Cb: Optional[npt.NDArray] = None
    U: Optional[npt.NDArray] = None
    V: Optional[npt.NDArray] = None
    P: Optional[npt.NDArray] = None
    free: Optional[List[int]] = field(default_factory=list)
    fixed: Optional[List[int]] = field(default_factory=list)
    ph: Optional[npt.NDArray] = None
    lb: Optional[npt.NDArray] = None
    ub: Optional[npt.NDArray] = None
    lb0: Optional[npt.NDArray] = None
    ub0: Optional[npt.NDArray] = None
    s: Optional[npt.NDArray] = None
    X: Optional[npt.NDArray] = None
    x0: Optional[npt.NDArray] = None
    y0: Optional[npt.NDArray] = None
    free_x: Optional[List[int]] = field(default_factory=list)
    free_y: Optional[List[int]] = field(default_factory=list)
    rol_x: Optional[List[int]] = field(default_factory=list)
    rol_y: Optional[List[int]] = field(default_factory=list)
    Citx: Optional[npt.NDArray] = None
    City: Optional[npt.NDArray] = None
    Cbtx: Optional[npt.NDArray] = None
    Cbty: Optional[npt.NDArray] = None
    xlimits: Optional[npt.NDArray] = None
    ylimits: Optional[npt.NDArray] = None
    qmin: Optional[npt.NDArray] = None
    qmax: Optional[npt.NDArray] = None
    k_i: Optional[Dict[Any, int]] = field(default_factory=dict)
    uv_i: Optional[Dict[Any, int]] = field(default_factory=dict)
    i_uv: Optional[Dict[int, Any]] = field(default_factory=dict)
    ind: Optional[List[int]] = field(default_factory=list)
    k: Optional[int] = None
    dep: Optional[List[int]] = field(default_factory=list)
    B: Optional[npt.NDArray] = None
    variables: Optional[List[str]] = field(default_factory=list)
    features: Optional[List[str]] = field(default_factory=list)
    constraints: Optional[List[str]] = field(default_factory=list)
    d: Optional[Any] = None
    # shape: Optional["Shape"] = None
    thk: Optional[float] = None
    rho: Optional[float] = None
    min_lb: Optional[float] = 0.0
    envelope: Optional["Envelope"] = None

    ub_lb_update: Optional[Callable] = None  # TODO: This needs to be taken care by the SurfaceModel
    b_update: Optional[Callable] = None  # TODO: This needs to be taken care by the SurfaceModel
    dub_dlb_update: Optional[Callable] = None  # TODO: This needs to be taken care by the SurfaceModel
    db_update: Optional[Callable] = None  # TODO: This needs to be taken care by the SurfaceModel
    Pmatrix: Optional[npt.NDArray] = None

    # Additional attributes for fixed problems
    Edinv: Optional[npt.NDArray] = None
    Ed: Optional[npt.NDArray] = None
    Ei: Optional[npt.NDArray] = None
    d0: Optional[npt.NDArray] = None
    time_inds: Optional[float] = None

    # Optimization functions (callables)
    fobj: Optional[Callable] = None
    fgrad: Optional[Callable] = None
    fconstr: Optional[Callable] = None
    fjac: Optional[Callable] = None

    # Optimization variables and initial values
    x0: Optional[npt.NDArray] = None
    bounds: Optional[List] = None
    f0: Optional[float] = None
    g0: Optional[npt.NDArray] = None

    @classmethod
    def from_formdiagram(cls, form: FormDiagram) -> "Problem":
        """Create a Problem from a FormDiagram.

        This is the primary constructor for creating a Problem instance. It initializes
        all matrices and vectors needed for the optimization from the form diagram.

        Parameters
        ----------
        form : :class:`~compas_tna.diagrams.FormDiagram`
            The FormDiagram containing the network geometry, loads, and constraints.

        Returns
        -------
        :class:`Problem`
            A fully initialized Problem instance with all necessary matrices.

        Examples
        --------
        >>> problem = Problem.from_formdiagram(form)
        >>> problem.m  # number of edges
        >>> problem.n  # number of vertices

        """
        # Mapping
        k_i = form.vertex_index()
        uv_i = form.uv_index()
        i_uv = form.index_uv()

        # Vertices and edges
        n = form.number_of_vertices()
        m = len(list(form.edges_where({"_is_edge": True})))
        fixed = [k_i[key] for key in form.supports()]
        nb = len(fixed)
        edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({"_is_edge": True})]
        free = list(set(range(n)) - set(fixed))
        ni = len(free)

        q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)
        qmax = array([form.edge_attribute((u, v), "qmax") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)
        qmin = array([form.edge_attribute((u, v), "qmin") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)

        # Co-ordinates, loads and constraints
        xyz = zeros((n, 3))
        x = zeros((n, 1))
        y = zeros((n, 1))
        z = zeros((n, 1))
        s = zeros((n, 1))  # target height
        px = zeros((n, 1))
        py = zeros((n, 1))
        pz = zeros((n, 1))
        lb = zeros((n, 1))
        ub = zeros((n, 1))
        xlimits = zeros((n, 2))
        ylimits = zeros((n, 2))

        for key, vertex in form.vertex.items():
            i = k_i[key]
            xyz[i, :] = form.vertex_coordinates(key)
            x[i] = vertex.get("x")
            y[i] = vertex.get("y")
            z[i] = vertex.get("z")
            px[i] = vertex.get("px", 0)
            py[i] = vertex.get("py", 0)
            pz[i] = vertex.get("pz", 0)
            s[i] = vertex.get("target", 0)  # used for bestfit
            if abs(s[i]) < 1e-6:
                s[i] = 0.0
            xlimits[i, 0] = vertex.get("xmin", None)
            xlimits[i, 1] = vertex.get("xmax", None)
            ylimits[i, 0] = vertex.get("ymin", None)
            ylimits[i, 1] = vertex.get("ymax", None)
            lb[i] = vertex.get("lb", None)
            ub[i] = vertex.get("ub", None)

        # Partial supports (rollers)
        rol_x = []
        rol_y = []

        for key in form.vertices_where({"rol_x": True}):
            rol_x.append(k_i[key])
        for key in form.vertices_where({"rol_y": True}):
            rol_y.append(k_i[key])

        free_x = list(set(free) - set(rol_x))
        free_y = list(set(free) - set(rol_y))

        # C and E matrices
        C = connectivity_matrix(edges, "csr")
        Ci = C[:, free]
        Cb = C[:, fixed]
        Cbtx = C[:, rol_x].transpose()
        Cbty = C[:, rol_y].transpose()
        Ct = C.transpose()
        Cit = Ci.transpose()
        Citx = Ct[free_x, :]
        City = Ct[free_y, :]
        uvw = C.dot(xyz)
        U = diags(uvw[:, 0].flatten())
        V = diags(uvw[:, 1].flatten())
        E = svstack((Citx.dot(U), City.dot(V))).toarray()

        # Settings of a free to move problem 'all q are variable'
        ind = list(range(m))
        dep = []
        k = m
        B = identity(m)
        d = zeros((m, 1))

        # Permutation matrix to fixed and free keys
        # with this matrix z = [Pmatrix]z will have the first ni nodes as free
        # and the last nb nodes as fixed
        Pfree = zeros((n, ni))
        Pfixed = zeros((n, nb))
        Pfree[free, :] = identity(ni)
        Pfixed[fixed, :] = identity(nb)
        Pmatrix = hstack([Pfree, Pfixed])

        # Create Problem instance
        return cls(
            q=q,
            m=m,
            n=n,
            ni=ni,
            nb=nb,
            E=E,
            C=C,
            Ct=Ct,
            Ci=Ci,
            Cit=Cit,
            Cb=Cb,
            U=U,
            V=V,
            P=hstack([px, py, pz]),
            free=free,
            fixed=fixed,
            ph=vstack([px[free_x], py[free_y]]),
            lb=lb,
            ub=ub,
            lb0=lb,
            ub0=ub,
            s=s,
            X=hstack([x, y, z]),
            x0=x,
            y0=y,
            free_x=free_x,
            free_y=free_y,
            rol_x=rol_x,
            rol_y=rol_y,
            Citx=Citx,
            City=City,
            Cbtx=Cbtx,
            Cbty=Cbty,
            xlimits=xlimits,
            ylimits=ylimits,
            qmin=qmin,
            qmax=qmax,
            k_i=k_i,
            uv_i=uv_i,
            i_uv=i_uv,
            ind=ind,
            k=k,
            dep=dep,
            B=B,
            d=d,
            Pmatrix=Pmatrix,
        )

    def _to_specialized(self, specialized_class, **updates):
        """Convert this Problem to a specialized subclass.

        This helper method creates a new instance of a specialized Problem class
        by copying all attributes from this base Problem and updating specific fields.

        Parameters
        ----------
        specialized_class : type
            The specialized class (FixedProblem, SymmetricProblem, etc.)
        **updates
            Fields to add/override for the specialized problem.

        Returns
        -------
        Problem subclass
            New instance of the specialized class with all base attributes plus updates.

        Examples
        --------
        >>> base = Problem.from_formdiagram(form)
        >>> fixed = base._to_specialized(FixedProblem, ind=ind, dep=dep, B=B, d=d)

        """
        from dataclasses import asdict

        # Get all fields from this problem instance
        data = asdict(self)

        # Update with specialized fields
        data.update(updates)

        # Create new instance of specialized class
        return specialized_class(**data)


# =============================================================================
# Specialized Problem Classes
# =============================================================================


@dataclass
class FixedProblem(Problem):
    """Problem with fixed diagram where independent edges are computed.

    This specialized problem class represents a structural analysis problem where
    the form diagram is fixed in plan, requiring the computation of independent
    edges for the force density method.

    The class inherits all attributes from Problem and automatically computes
    the independent edges and associated matrices during initialization.

    """

    @classmethod
    def from_formdiagram(
        cls,
        form: FormDiagram,
        method: str = "QR",
        printout: bool = False,
        tol: Optional[float] = None,
    ) -> "FixedProblem":
        """Create a FixedProblem from a FormDiagram with independent edges computed.

        Parameters
        ----------
        form : :class:`~compas_tna.diagrams.FormDiagram`
            The FormDiagram containing the network geometry, loads, and constraints.
        method : str, optional
            Method to find independent edges, by default 'QR'.
            Options: 'QR', 'SVD', etc.
        printout : bool, optional
            Whether to print progress information, by default False.
        tol : float, optional
            Tolerance for the singular values when finding independents, by default None.

        Returns
        -------
        :class:`FixedProblem`
            A fully initialized FixedProblem with independent edges computed.

        Examples
        --------
        >>> form = FormDiagram.create_cross(x_span=(0, 10), y_span=(0, 10), n=5)
        >>> problem = FixedProblem.from_formdiagram(form, method="QR")
        >>> problem.k  # number of independent edges (< m)

        """
        # Create base problem
        base_problem = Problem.from_formdiagram(form)

        # Find independent edges
        ind = []
        start_time = time.time()

        ind_edges = form.edges_where({"_is_ind": True})

        if len(list(ind_edges)) > 0:
            ind = [base_problem.uv_i[edge] for edge in ind_edges]
        else:
            ind = find_independents(base_problem.E, method=method, tol=tol)

        k = len(ind)
        dep = list(set(range(base_problem.m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print("Reduced problem to {0} force variables with ind. edges".format(k))
            print("Elapsed Time: {0:.1f} sec".format(elapsed_time))

        # Update form diagram with independent edge information
        points = []
        for edge in form.edges_where({"_is_edge": True}):
            if base_problem.uv_i[edge] in ind:
                form.edge_attribute(edge, "is_ind", True)
                points.append(Point(*(form.edge_midpoint(edge)[:2] + [0])))
            else:
                form.edge_attribute(edge, "is_ind", False)
        form.attributes["indset"] = points

        # Compute matrices for fixed problem
        rcond = 1e-17
        if tol:
            rcond = tol
        Ed = base_problem.E[:, dep]
        Edinv = -csr_matrix(pinv(base_problem.E[:, dep], rcond=rcond))
        Ei = csr_matrix(base_problem.E[:, ind])
        B = zeros((base_problem.m, k))
        B[dep] = Edinv.dot(Ei).toarray()
        B[ind] = identity(k)

        d = zeros((base_problem.m, 1))
        d[dep] = -Edinv.dot(base_problem.ph)  # q = Bqi + d | d = Ed(-1)*ph

        if any(base_problem.ph):
            check_hor = check_horizontal_loads(base_problem.E, base_problem.ph)
            if check_hor:
                if printout:
                    print("Horizontal Loads can be taken!")
            else:
                print("Warning: Horizontal Loads are not suitable for this FD!")

        # Create FixedProblem from base with updated attributes
        return base_problem._to_specialized(
            cls,
            ind=ind,
            k=k,
            dep=dep,
            B=B,
            d=d,
            d0=d,
            Edinv=Edinv,
            Ed=Ed,
            Ei=Ei,
            time_inds=elapsed_time,
        )


@dataclass
class SymmetricProblem(Problem):
    """Problem with symmetry constraints applied.

    This specialized problem class represents a structural analysis problem where
    the form diagram has symmetry, allowing reduction of variables by exploiting
    symmetric patterns.

    """

    @classmethod
    def from_formdiagram(
        cls,
        form: FormDiagram,
        list_axis_symmetry: Optional[list] = None,
        center=None,
        correct_loads: bool = True,
        printout: bool = False,
    ) -> "SymmetricProblem":
        """Create a SymmetricProblem from a FormDiagram with symmetry constraints.

        Parameters
        ----------
        form : :class:`~compas_tna.diagrams.FormDiagram`
            The FormDiagram containing the network geometry, loads, and constraints.
        list_axis_symmetry : list, optional
            List of the axis of symmetry to consider, by default None.
            If None, symmetry is inferred from the diagram type.
        center : list, optional
            The center of the pattern if it's a circular pattern, by default None.
        correct_loads : bool, optional
            Whether to update applied loads regarding the symmetry, by default True.
        printout : bool, optional
            Whether to print progress information, by default False.

        Returns
        -------
        :class:`SymmetricProblem`
            A fully initialized SymmetricProblem with symmetry applied.

        Examples
        --------
        >>> form = FormDiagram.create_circular_radial(...)
        >>> problem = SymmetricProblem.from_formdiagram(form)
        >>> problem.k  # number of independent edges (reduced by symmetry)

        """
        # Create base problem
        base_problem = Problem.from_formdiagram(form)

        start_time = time.time()

        # Apply symmetry to form
        apply_sym_to_form(form, list_axis_symmetry, center, correct_loads)

        # Build symmetry transformation
        Esym = build_symmetry_transformation(form, printout=False)
        mapsym = form.build_symmetry_map()
        ind = sorted(list(mapsym.values()))

        k = len(ind)
        dep = list(set(range(base_problem.m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print("Reduced problem to {0} force variables by SYM".format(k))
            print("Elapsed Time: {0:.1f} sec".format(elapsed_time))

        # Mark independent edges on form
        for u, v in form.edges_where({"_is_edge": True}):
            form.edge_attribute((u, v), "is_ind", True if base_problem.uv_i[(u, v)] in ind else False)

        B = Esym

        # Create SymmetricProblem from base with updated attributes
        return base_problem._to_specialized(
            cls,
            ind=ind,
            k=k,
            dep=dep,
            B=B,
        )


@dataclass
class FixedSymmetricProblem(Problem):
    """Problem with both fixed diagram and symmetry constraints.

    This specialized problem class combines the features of FixedProblem and
    SymmetricProblem, applying both independent edge computation and symmetry
    reduction.

    """

    @classmethod
    def from_formdiagram(
        cls,
        form: FormDiagram,
        method: str = "QR",
        list_axis_symmetry: Optional[list] = None,
        center=None,
        correct_loads: bool = True,
        printout: bool = False,
        tol: Optional[float] = None,
    ) -> "FixedSymmetricProblem":
        """Create a FixedSymmetricProblem from a FormDiagram.

        Parameters
        ----------
        form : :class:`~compas_tna.diagrams.FormDiagram`
            The FormDiagram containing the network geometry, loads, and constraints.
        method : str, optional
            Method to find independent edges, by default 'QR'.
        list_axis_symmetry : list, optional
            List of the axis of symmetry to consider, by default None.
        center : list, optional
            The center of the pattern if it's a circular pattern, by default None.
        correct_loads : bool, optional
            Whether to update applied loads regarding the symmetry, by default True.
        printout : bool, optional
            Whether to print progress information, by default False.
        tol : float, optional
            Tolerance for the singular values, by default None.

        Returns
        -------
        :class:`FixedSymmetricProblem`
            A fully initialized problem with both adaptations applied.

        Examples
        --------
        >>> form = FormDiagram.create_circular_radial(...)
        >>> problem = FixedSymmetricProblem.from_formdiagram(form, method="QR")
        >>> problem.k  # number of independent edges (reduced by both)

        """
        # First, create FixedProblem (this handles independent edge computation)
        fixed_problem = FixedProblem.from_formdiagram(form, method=method, printout=False, tol=tol)

        start_time = time.time()

        # Then, apply symmetry reduction on top of fixed problem
        apply_sym_to_form(form, list_axis_symmetry, center, correct_loads)

        i_uv = fixed_problem.i_uv
        uv_i = fixed_problem.uv_i
        unique_sym = []
        ind_reduc = []

        for index in fixed_problem.ind:
            u, v = i_uv[index]
            index_sym = form.edge_attribute((u, v), "sym_key")
            if index_sym not in unique_sym:
                unique_sym.append(index_sym)
                ind_reduc.append(index)

        k_sym = len(ind_reduc)
        Bsym = zeros((fixed_problem.k, k_sym))

        j = 0
        for key in unique_sym:
            edges_sym = [uv_i[(u, v)] for u, v in form.edges_where({"sym_key": key})]
            col = zeros((fixed_problem.k, 1))
            for i, ind_i in enumerate(fixed_problem.ind):
                if ind_i in edges_sym:
                    col[i] = 1.0
            Bsym[:, [j]] = col
            j += 1

        # Combine fixed and symmetry transformations
        B_combined = fixed_problem.B.dot(Bsym)
        ind_final = ind_reduc
        k_final = k_sym
        dep_final = list(set(range(fixed_problem.m)) - set(ind_final))

        elapsed_time = time.time() - start_time

        if printout:
            print("Reduced problem to {0} force variables".format(k_final))
            print("Elapsed Time: {0:.1f} sec".format(elapsed_time))

        # Mark final independent edges
        for u, v in form.edges_where({"_is_edge": True}):
            form.edge_attribute((u, v), "is_ind", True if fixed_problem.uv_i[(u, v)] in ind_final else False)

        # Create FixedSymmetricProblem from FixedProblem with symmetry updates
        return fixed_problem._to_specialized(
            cls,
            ind=ind_final,
            k=k_final,
            dep=dep_final,
            B=B_combined,
            time_inds=fixed_problem.time_inds + elapsed_time,
        )


# =============================================================================
# Helpers
# =============================================================================


def apply_sym_to_form(form: FormDiagram, list_axis_symmetry=None, center=None, correct_loads=True) -> None:
    """Apply symmetry to the form diagram.

    Parameters
    ----------
    problem : :class:`~compas_tno.problems.Problem`
        The problem with matrices for calculation
    form : :class:`~compas_tna.diagrams.FormDiagram`
        The form diagram to analyse
    list_axis_symmetry : [list], optional
        List of the axis of symmetry to consider, by default None
    center : [list], optional
        The center of the pattern if it's a circula pattern, by default None
    correct_loads : bool, optional
        If update should be done in the applied loads regarding the symmetry, by default True

    """
    if not list_axis_symmetry:
        data = form.parameters
        if data["type"] in ["radial_fd", "radial_spaced_fd", "spiral_fd"]:
            apply_radial_symmetry(form, center=center, correct_loads=correct_loads)
        elif data["type"] in ["cross_fd", "fan_fd", "cross_diagonal", "cross_with_diagonal", "ortho_fd"]:
            list_axis_symmetry = find_sym_axis_in_rect_patterns(data)
            print("Axis of Symmetry identified:", list_axis_symmetry)
        else:
            print("Symmetry applied to a unknown form diagram type")
            raise NotImplementedError

    if list_axis_symmetry:
        apply_symmetry_from_axis(form, list_axis_symmetry=list_axis_symmetry, correct_loads=correct_loads)


def check_bad_independents(problem: Problem, form: FormDiagram, printout=False, eps=10e6) -> None:
    keep_inds = []

    for i_ind in range(len(problem.ind)):
        column_i = abs(problem.B[:, i_ind]).flatten()
        max_col = max(column_i)
        if max_col < eps:
            keep_inds.append(i_ind)

    if len(keep_inds) == problem.m:
        return

    print("original inds:", len(problem.ind), problem.ind)
    print("Updating the relevant independent")
    inds = []
    for i_ind in keep_inds:
        inds.append(problem.ind[i_ind])

    q0 = zeros((problem.m, 1))
    form.edges_attribute("is_ind", False)
    for i, edge in enumerate(form.edges_where({"_is_edge": True})):
        if i in inds:
            form.edge_attribute(edge, "is_ind", True)
        elif i in problem.ind:
            q0[i] = form.edge_attribute(edge, "q")

    print(q0)

    problem.ind = inds
    problem.dep = list(set(range(problem.m)) - set(inds))
    problem.B = problem.B[:, keep_inds]
    problem.d += q0
    problem.k = len(inds)
    print("final inds", len(problem.ind), problem.ind)


# def plot_svds(M, tol=None):
#     """Plot a diagram of SVD.

#     Parameters
#     ----------
#     M : :class:`~compas_tno.problems.Problem`
#         The problem with matrices for calculation
#     tol : float (optional)
#         The tolerance to find the independents
#     """

#     import matplotlib.pyplot as plt
#     from numpy.linalg import matrix_rank

#     k = len(M.ind)
#     n, m = M.E.shape
#     # mn_min = min(n, m)

#     _, s, _ = svd(asarray(M.E))
#     # _, s, _ = svds(M.E, k=min(M.E.shape), solver='propack')
#     print("max/min singular vectors E", max(s), min(s), len(s))
#     print("Shape E: {} | #ind: {} | rank : {}:".format(M.E.shape, k, matrix_rank(M.E, tol=tol)))

#     # mn = max(m - n, n - m)
#     mn = max(m - n, 0)
#     zs = k - mn
#     print("# null-SVD:", zs)

#     # if zs + 3 > len(s):
#     #     print('Last {} SVs: {}'.format(len(s[-(zs + 3):]), s[-(zs + 3):]))
#     # print('ALL SVS:', s)

#     check = check_independents(M)
#     print("Check independents:", check)
#     print("Shape/ Rank Ed:", M.Ed.shape, matrix_rank(M.Ed))

#     if zs > 0:
#         first_zero = s[len(s) - zs]
#         last_non_zero = s[len(s) - zs - 1]

#         x_zs = list(range(len(s)))[-zs:]
#         y_zs = s[-zs:]

#         print("Last Non Zero:", last_non_zero)
#         print("First Zero SV:", first_zero)

#         # len_zero = len(s[s<1.0])
#         lin_x = len(s) - zs

#         porcentage_key = (last_non_zero - first_zero) / last_non_zero
#         print("percentage key is:", porcentage_key)

#         _, ax = plt.subplots()
#         ax.plot(s)
#         ax.scatter(x_zs, y_zs)
#         ax.plot([lin_x, lin_x], [0, 1.0], color="black")
#         ax.plot()
#         plt.show()

#     else:
#         print("EQ. Matrix is FULL rank: (k = m-n)")
#         _, ax = plt.subplots()
#         # ax.set_ylim(2.7, 3.0)
#         ax.plot(s)
#         plt.show()

#     return
