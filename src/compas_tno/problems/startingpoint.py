from compas_tno.algorithms import compute_reactions
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.solvers import run_loadpath_from_form_CVXPY


def startingpoint_sag(form, boundary_force=50.0, **kwargs):
    """Initialize the equilibrium in a form diagram with applied loads using sag approach

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram. Loads and support must already have been assigned
    boundary_force : float, optional
        Force density in the edges on the boundary.
        The default value is ``50.0``.
    """
    form.edges_attribute("q", min(boundary_force, -1 * boundary_force), list(form.edges_on_boundary()))
    startingpoint_fdm(form)
    for key in form.vertices():
        form.vertex_attribute(key, "z", 0.0)
    return form


def startingpoint_loadpath(form, problem=None, find_inds=False, solver_convex="CLARABEL", printout=False, **kwargs):
    """Built-in function to optimise the loadpath considering diagram fixed projection.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to compute the loadpath
    problem : :class:`~compas_tno.problems.Problem`, optional
        The problem classs with the matrices relevant of the problem, by default None
    find_inds : bool, optional
        If independents need to be found before the loadpath computation, by default False
    solver_convex : str, optional
        Solver to use, by default CLARABEL. Options are "CLARABEL", "MOSEK" or "CVXOPT".
        Note: "MOSEK" and "CVXOPT" are not available in the default installation of TNO.
    printout : bool, optional
        If prints about the optimisation data should appear in the screen, by default False

    Returns
    -------
    problem : Problem
        The class with the main matrices of the problem
    """

    problem = run_loadpath_from_form_CVXPY(
        form,
        problem=problem,
        find_inds=find_inds,
        solver_convex=solver_convex,
        printout=printout,
    )

    return problem


def startingpoint_tna(form, plot=False, **kwargs):
    """Initialize the equilibrium in a form diagram with applied loads using TNA interative solver procedure (form and force diagrams are parallel)

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram. Loads and support must already have been assigned
    plot : bool, optional
        Plots of the intermediare and final force diagrams to follow the process, by default False
    """

    # form_update_with_parallelisation(form, plot=plot)

    # compute_reactions(form)

    raise NotImplementedError("Starting point not implemented for TNA")


# Use the appropiate functions at TNA here
def startingpoint_fdm(form, **kwargs):
    """Initialize the equilibrium in a form diagram with applied loads using FDM approach for the q's stored in the form

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram. Loads, supports and force densities must already have been assigned
    """

    equilibrium_fdm(form)

    compute_reactions(form)
