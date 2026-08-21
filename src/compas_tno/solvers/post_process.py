from typing import TYPE_CHECKING
from typing import Optional

if TYPE_CHECKING:
    from compas_tna.diagrams import FormDiagram
    from compas_tno.analysis import Analysis
    from compas_tno.problems import Problem
    from compas_tno.solvers.result import SolverResult

from compas_tno.algorithms import compute_reactions
from compas_tno.algorithms import q_from_variables
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import save_geometry_at_iterations


def apply_solution_to_form(
    form: "FormDiagram",
    problem: "Problem",
    result: "SolverResult",
    settings: Optional[dict] = None,
) -> "Problem":
    """Apply optimization result to form diagram and problem.

    This is the clean, decoupled version of post-processing that uses
    the SolverResult directly without depending on Analysis or Optimiser.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to update.
    problem : Problem
        The problem structure with matrices and functions.
    result : SolverResult
        The optimization result from solver.
    settings : dict, optional
        Settings for post-processing:
        - 'printout' : bool - Print summary (default: True)
        - 'summary' : bool - Print detailed summary (default: False)
        - 'save_iterations' : bool - Save iteration data (default: False)
        - 'save_force_diagram' : bool - Save force diagram (default: True)

    Returns
    -------
    Problem
        The updated problem with optimal solution applied.

    Examples
    --------
    >>> result = solver.solve(problem)
    >>> problem = apply_solution_to_form(form, problem, result)
    >>> # form and problem are now updated with optimal solution

    """
    settings = settings or {}
    printout = settings.get("printout", True)
    summary = settings.get("summary", False)

    # Extract results from SolverResult
    xopt = result.xopt
    fopt = result.fopt
    message = result.message

    # Get constraint function from problem
    fconstr = problem.fconstr
    thk = problem.thk
    i_uv = form.index_uv()

    # Extract variables from optimal solution
    check = problem.k
    qid = xopt[: problem.k].reshape(-1, 1)

    if "xyb" in problem.variables:
        xyb = xopt[check : check + 2 * problem.nb]
        check = check + 2 * problem.nb
        problem.X[problem.fixed, :2] = xyb.reshape(-1, 2, order="F")
    if "zb" in problem.variables:
        zb = xopt[check : check + problem.nb]
        check = check + problem.nb
        problem.X[problem.fixed, [2]] = zb.flatten()
    if "t" in problem.variables:
        thk = xopt[-1]
    if "n" in problem.variables:
        n = xopt[-1]
    if "lambdh" in problem.variables:
        lambdh = xopt[check : check + 1]
        problem.P[:, [0]] = lambdh * problem.px0
        problem.P[:, [1]] = lambdh * problem.py0
        problem.d = lambdh * problem.d0
    if "lambdv" in problem.variables:
        lambdv = xopt[check : check + 1]
        problem.P[:, [2]] = lambdv * problem.pzv + problem.pz0
    if "tub" in problem.variables:
        tub = xopt[check : check + problem.n]
        check = check + problem.n
    if "tlb" in problem.variables:
        tlb = xopt[check : check + problem.n]
        check = check + problem.n
    if "tub_reac" in problem.variables:
        tub_reac = xopt[check : check + 2 * problem.nb]
        check = check + 2 * problem.nb

    # Update problem with optimal force densities
    problem.q = q_from_variables(qid, problem.B, problem.d)

    # Compute final constraint values
    if fconstr:
        g_final = fconstr(xopt, problem)
    else:
        g_final = None

    # Update geometry
    problem.X[problem.free] = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)

    # Update form diagram vertices
    for i, key in enumerate(form.vertices()):
        form.vertex_attribute(key, "x", problem.X[i, 0])
        form.vertex_attribute(key, "y", problem.X[i, 1])
        form.vertex_attribute(key, "z", problem.X[i, 2])
        form.vertex_attribute(key, "px", problem.P[i, 0])
        form.vertex_attribute(key, "py", problem.P[i, 1])
        form.vertex_attribute(key, "pz", problem.P[i, 2])

    # Update form diagram edges
    for c, qi in enumerate(list(problem.q.ravel())):
        u, v = i_uv[c]
        li = form.edge_length((u, v))
        form.edge_attribute((u, v), "q", float(qi))
        form.edge_attribute((u, v), "_f", float(qi * li))

    # Compute reactions
    compute_reactions(form)

    # Handle thickness updates
    if "t" in problem.variables:
        if problem.envelope.is_parametric:
            problem.envelope.thickness = thk
            problem.envelope.update_envelope()
            problem.envelope.apply_bounds_to_formdiagram(form)
            if printout:
                print(f"Envelope updated for new minimum thickness: {thk}")
        else:
            pass
        # # TODO: Now that only the model is passed, if the optimisation is minimum thickness, we need to update the model based on the template
        # if shape.parameters["type"] == "general":
        #     if thickness_type == "constant":
        #         form.attributes["thk"] = thk
        #         shape.parameters["thk"] = thk
        #         shape.intrados = shape.middle.offset_mesh(n=thk / 2, direction="down")
        #         shape.extrados = shape.middle.offset_mesh(n=thk / 2, direction="up")
        #         form.apply_envelope_from_shape(shape)

        #     elif thickness_type == "variable":
        #         t0 = shape.parameters["thk"]
        #         thk = t0 * thk  # Consider that the thk for general shapes is a percentage of the thickness
        #         form.attributes["thk"] = thk
        #         shape.parameters["thk"] = thk
        #         if printout:
        #             print("Optimum Value corresponds to a thickness of:", thk)
        #         shape.extrados, shape.intrados = shape.middle.offset_up_and_down(n=fopt)
        #         form.apply_envelope_from_shape(shape)

        #     elif thickness_type == "intrados":
        #         form.attributes["thk"] = thk
        #         shape.parameters["thk"] = thk
        #         shape.middle = shape.intrados.offset_mesh(n=thk / 2, direction="up")
        #         shape.extrados = shape.intrados.offset_mesh(n=thk, direction="up")
        #         form.apply_envelope_from_shape(shape)
        # else:
        #     form.attributes["thk"] = thk
        #     shape.parameters["thk"] = thk
        #     shape = Shape.from_library(shape.parameters)
        #     form.apply_envelope_from_shape(shape)  # Check if this is ok for adapted pattern
        #     form.apply_bounds_reactions(shape)

        #     i = 0
        #     for key in form.vertices():  # this resolve the problem due to the adapted pattern
        #         form.vertex_attribute(key, "ub", float(problem.ub[i]))
        #         form.vertex_attribute(key, "lb", float(problem.lb[i]))
        #         i += 1

    # if "update-envelope" in features:
    #     form.attributes["thk"] = thk
    #     shape.parameters["thk"] = thk
    #     shape = Shape.from_library(shape.parameters)
    #     form.apply_envelope_from_shape(shape)

    # if 's' in problem.variables:
    #     s = -1 * fopt
    #     for key in form.vertices():
    #         ub = form.vertex_attribute(key, 'ub')
    #         lb = form.vertex_attribute(key, 'lb')
    #         form.vertex_attribute(key, 'ub', ub - s * (ub - lb))
    #         form.vertex_attribute(key, 'lb', lb + s * (ub - lb))

    # Handle variable-specific updates
    if "n" in problem.variables:
        if printout:
            print(f"Value of N: {n}")

    if "tub" in problem.variables:
        for i, key in enumerate(form.vertices()):
            zub = form.vertex_attribute(key, "ub")
            form.vertex_attribute(key, "tub", tub[i])
            form.vertex_attribute(key, "ub", zub + tub[i])

    if "tlb" in problem.variables:
        for i, key in enumerate(form.vertices()):
            zlb = form.vertex_attribute(key, "lb")
            form.vertex_attribute(key, "tlb", tlb[i])
            form.vertex_attribute(key, "lb", zlb - tlb[i])

    if "tub_reac" in problem.variables:
        for i, key in enumerate(form.vertices_where({"is_support": True})):
            form.vertex_attribute(key, "tub_reac", [tub_reac[i], tub_reac[i + problem.nb]])

    if "lambdv" in problem.variables:
        added_load = lambdv * problem.pzv
        for i, key in enumerate(form.vertices()):
            p_added = added_load[i]
            form.vertex_attribute(key, "pzext", float(p_added))

    if "lambdh" in problem.variables:
        form.attributes["lambdh"] = lambdh

    # Print summary
    if printout or summary:
        print("\n" + "-" * 50)
        print("TNO v.2.0")
        print(f"Solution  : {message}")
        print(f"Status    : {'SUCCESS' if result.success else 'FAILED'}")
        try:
            print(f"q range   : {min(problem.q):.3f} : {max(problem.q):.3f}")
        except BaseException:
            print(f"q range   : {min(problem.q.flatten()):.3f} : {max(problem.q.flatten()):.3f}")
        print(f"zb range  : {min(problem.X[problem.fixed, [2]]):.3f} : {max(problem.X[problem.fixed, [2]]):.3f}")
        if g_final is not None:
            print(f"constr    : {min(g_final):.3f} : {max(g_final):.3f}")
        print(f"fopt      : {fopt:.3f}")
        if result.niter is not None:
            print(f"iterations: {result.niter}")
        print(f"time      : {result.time:.3f} sec")
        print("-" * 50 + "\n")

    return problem


def post_process_nlopt(analysis: "Analysis"):
    """Post processing of the optimisation.

    This function is kept for backward compatibility with the Analysis workflow.
    It wraps the new apply_solution_to_form() function.

    Parameters
    ----------
    analysis : Analysis
        The Analysis object

    Returns
    -------
    analysis : Analysis
        The Analysis object updated

    Notes
    -----
    For new code, prefer using::

        result = solver.solve(problem)
        problem = apply_solution_to_form(form, problem, result, settings)

    """
    form = analysis.formdiagram
    optimiser = analysis.optimiser
    problem = optimiser.problem

    # Get SolverResult from analysis
    result = analysis.result

    # Call the new clean function
    problem = apply_solution_to_form(
        form=form,
        problem=problem,
        result=result,
        settings=optimiser.settings,
    )

    # Handle save_iterations (needs optimiser for legacy reasons)
    if optimiser.settings.get("save_iterations", False):
        show_force_diagram = optimiser.settings.get("save_force_diagram", True)
        file_Xform, file_Xforce = save_geometry_at_iterations(form, optimiser, force=show_force_diagram)
        analysis.optimiser.Xform = file_Xform
        analysis.optimiser.Xforce = file_Xforce

    return analysis
