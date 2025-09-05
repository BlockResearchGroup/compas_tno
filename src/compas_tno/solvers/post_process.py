from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from compas_tna.diagrams import FormDiagram
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.problems import Problem

from compas_tno.algorithms import compute_reactions
from compas_tno.algorithms import q_from_variables
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import save_geometry_at_iterations


def post_process_nlopt(analysis: "Analysis"):
    """Post processing of the optimisation.

    Parameters
    ----------
    analysis : Analysis
        The Analysis object

    Returns
    -------
    analysis : Analysis
        The Analysis object updated
    """

    form: "FormDiagram" = analysis.formdiagram
    optimiser: "Optimiser" = analysis.optimiser
    # envelope: "Envelope" = analysis.envelope

    problem: "Problem" = optimiser.problem
    summary = optimiser.settings.get("summary", False)
    printout = optimiser.settings.get("printout", True)
    # thickness_type = optimiser.settings.get("thickness_type", "constant")
    # features = optimiser.settings.get("features", [])
    save_iterations = optimiser.settings.get("save_iterations", False)
    show_force_diagram = optimiser.settings.get("save_force_diagram", True)

    fconstr = optimiser.fconstr
    xopt = optimiser.xopt
    fopt = optimiser.fopt
    message = optimiser.message
    thk = problem.thk

    i_uv = form.index_uv()

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

    problem.q = q_from_variables(qid, problem.B, problem.d)

    # ADD 'lambdv' to post-process

    # if 's' in variables:
    #     s = xopt[-1]

    g_final = fconstr(xopt, problem)
    problem.X[problem.free] = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)

    # if printout:
    #     print('post-processing min, max z:', min(problem.X[:, 2].flatten()), max(problem.X[:, 2].flatten()))

    i = 0
    for key in form.vertices():
        form.vertex_attribute(key, "x", problem.X[i, 0])
        form.vertex_attribute(key, "y", problem.X[i, 1])
        form.vertex_attribute(key, "z", problem.X[i, 2])
        form.vertex_attribute(key, "px", problem.P[i, 0])
        form.vertex_attribute(key, "py", problem.P[i, 1])
        form.vertex_attribute(key, "pz", problem.P[i, 2])
        i = i + 1

    for c, qi in enumerate(list(problem.q.ravel())):
        u, v = i_uv[c]
        li = form.edge_length((u, v))
        form.edge_attribute((u, v), "q", float(qi))
        form.edge_attribute((u, v), "f", float(qi * li))

    # form.attributes["loadpath"] = form.loadpath()
    compute_reactions(form)

    if "t" in problem.variables:
        if problem.envelope.is_parametric:
            problem.envelope.thickness = thk
            problem.envelope.update_envelope()
            problem.envelope.apply_bounds_to_formdiagram(form)
            print("Envelope updated for new minimum thickness: {}".format(thk))
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

    if "n" in problem.variables:
        print("Value of N:", n)
        n = -1 * fopt
        pass
        # shape.intrados = shape.intrados.offset_mesh(n=n, direction="up")
        # shape.extrados = shape.extrados.offset_mesh(n=n, direction="down")
        # form.apply_envelope_from_shape(shape)

    if "tub" in problem.variables:
        for i, key in enumerate(form.vertices()):
            zub = form.vertex_attribute(key, "ub")
            form.vertex_attribute(key, "tub", tub[i])
            form.vertex_attribute(key, "ub", zub + tub[i])

    if "tlb" in problem.variables:
        for i, key in enumerate(form.vertices()):
            zub = form.vertex_attribute(key, "lb")
            form.vertex_attribute(key, "tlb", tlb[i])
            form.vertex_attribute(key, "lb", zub - tlb[i])

    if "tub_reac" in problem.variables:
        for i, key in enumerate(form.vertices_where({"is_support": True})):
            print(i, key, tub_reac)
            form.vertex_attribute(key, "tub_reac", [tub_reac[i], tub_reac[i + problem.nb]])

    if "lambdv" in problem.variables:
        added_load = lambdv * problem.pzv
        for i, key in enumerate(form.vertices()):
            p_added = added_load[i]
            form.vertex_attribute(key, "pzext", float(p_added))

    if "lambdh" in problem.variables:
        form.attributes["lambdh"] = lambdh  # can be improved. It currently takes fopt, but if loads at start are 0.1*SWT, the obj is lambd/0.1

    # analysis.form = form
    # analysis.optimiser = optimiser
    # analysis.shape = shape

    if save_iterations:
        file_Xform, file_Xforce = save_geometry_at_iterations(form, optimiser, force=show_force_diagram)
        analysis.optimiser.Xform = file_Xform
        analysis.optimiser.Xforce = file_Xforce

    if printout or summary:
        print("\n" + "-" * 50)
        print("TNO v.2.0")
        print("Solution  :", message)
        try:
            print("q range : {0:.3f} : {1:.3f}".format(min(problem.q), max(problem.q)))
        except BaseException:
            print("q range : {0:.3f} : {1:.3f}".format(min(problem.q.flatten()), max(problem.q.flatten())))
        print("zb range  : {0:.3f} : {1:.3f}".format(min(problem.X[problem.fixed, [2]]), max(problem.X[problem.fixed, [2]])))
        print("constr    : {0:.3f} : {1:.3f}".format(min(g_final), max(g_final)))
        print("fopt      : {0:.3f}".format(fopt))
        print("-" * 50 + "\n")

    return analysis
