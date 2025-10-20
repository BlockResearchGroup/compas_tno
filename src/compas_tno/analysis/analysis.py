from typing import Annotated
from typing import Literal
from typing import Optional

import numpy.typing as npt

from compas.data import Data
from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import Envelope
from compas_tno.optimisers import Optimiser
from compas_tno.problems import set_up_base_problem
from compas_tno.problems import set_up_convex_optimisation
from compas_tno.problems import set_up_general_optimisation
from compas_tno.solvers import SolverResult
from compas_tno.solvers import post_process_nlopt
from compas_tno.solvers import run_convex_optimisation
from compas_tno.solvers import run_nlopt_ipopt
from compas_tno.solvers import run_nlopt_scipy


class Analysis(Data):
    """Class to organise the data involved in an optimisation.

    The ``Analysis`` class unites the :class:`~compas_tno.diagrams.FormDiagram`, the :class:`~compas_tno.shapes.Shape`
    and the :class:`~compas_tno.optimisers.Optimiser` to compute the optimisation.

    Examples
    --------
    >>> from compas_tno.analysis import Analysis
    >>> from compas_tna.shapes import Envelope
    >>> from compas_tna.diagrams import FormDiagram
    >>> arch = Envelope.create_arch()
    >>> form = FormDiagram.create_arch()
    >>> analysis = Analysis.create_minthrust_analysis(form, arch)

    """

    formdiagram: FormDiagram
    envelope: Envelope
    optimiser: Optimiser
    result: SolverResult
    settings: dict
    name: str

    def __init__(
        self,
        formdiagram: FormDiagram,
        envelope: Envelope,
        optimiser: Optimiser,
        settings: Optional[dict] = None,
        name: Optional[str] = None,
        result: Optional[SolverResult] = None,
    ):
        name = name or "Analysis"

        super().__init__(name=name)

        self.settings = settings or {}
        self.formdiagram = formdiagram
        self.envelope = envelope
        self.optimiser = optimiser
        self.result = result  # SolverResult from optimization

    def __str__(self):
        tpl = "<Analysis with parameters: {} >".format(self.data)
        return tpl

    @property
    def __data__(self) -> dict:
        data = {
            "formdiagram": self.formdiagram,
            "envelope": self.envelope,
            "optimiser": self.optimiser,
            "settings": self.settings,
        }
        return data

    @classmethod
    def create_minthk_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        solver: str = "SLSQP",
        derivatives: bool = True,
    ) -> "Analysis":
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_minthk_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            solver=solver,
            derivatives=derivatives,
        )

        if printout:
            print("-" * 20)
            print("Minimum thickness analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_bestfit_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        solver: str = "SLSQP",
        derivatives: bool = True,
    ) -> "Analysis":
        """Create a bestfit analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        derivatives : bool, optional
            Whether or not derivatives should be considered, by default False

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_bestfit_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            solver=solver,
            derivatives=derivatives,
        )

        if printout:
            print("-" * 20)
            print("Minimum thickness analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_minthrust_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        solver: str = "SLSQP",
    ) -> "Analysis":
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_minthrust_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            solver=solver,
        )

        if printout:
            print("-" * 20)
            print("Minimum thrust analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_maxthrust_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        solver: str = "SLSQP",
    ) -> "Analysis":
        """Create a maximum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_maxthrust_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            solver=solver,
        )

        if printout:
            print("-" * 20)
            print("Maximium thrust analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_max_load_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        plot: bool = False,
        horizontal: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        solver: str = "IPOPT",
        derivatives: bool = True,
        load_direction: Optional[Annotated[npt.NDArray, Literal["2n, 1"]]] = None,
        max_lambd: float = 1.0,
    ) -> "Analysis":
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        plot : bool, optional
            Whether or load applied is horizontal, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        if horizontal:
            optimiser = Optimiser.create_max_horload_optimiser(
                printout=printout,
                plot=plot,
                max_iter=max_iter,
                starting_point=starting_point,
                solver=solver,
                derivatives=derivatives,
                load_direction=load_direction,
                max_lambd=max_lambd,
            )

            if printout:
                print("-" * 20)
                print("Max horizontal load analysis created")
                # print(optimiser)

        else:
            optimiser = Optimiser.create_max_vertload_optimiser(
                printout=printout,
                plot=plot,
                max_iter=max_iter,
                starting_point=starting_point,
                solver=solver,
                derivatives=derivatives,
                load_direction=load_direction,
                max_lambd=max_lambd,
            )

            if printout:
                print("-" * 20)
                print("Max vertical load analysis created")
                # print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_compl_energy_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        solver: str = "IPOPT",
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        support_displacement: Optional[Annotated[npt.NDArray, Literal["nb, 3"]]] = None,
        Emethod: str = "simplified",
    ) -> "Analysis":
        """Create a complementary energy analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        support_displacement : array [nb x 3], optional
            Vector with the displacement applied to the supports, by default None
        Emethod : str, optional
            Whether or not internal deformations should be considered, by default 'simplified' which considers the material rigid

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_compl_energy_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            support_displacement=support_displacement,
            Emethod=Emethod,
            solver=solver,
        )

        if printout:
            print("-" * 20)
            print("Complementary energy created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_quad_compl_energy_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        printout: bool = False,
        solver: str = "IPOPT",
        plot: bool = False,
        max_iter: int = 500,
        starting_point: str = "loadpath",
        support_displacement: Optional[Annotated[npt.NDArray, Literal["nb, 3"]]] = None,
        Emethod: str = "simplified",
    ):
        """Create a complementary energy analysis including a quadratic term from the elements of the problem (form and shape)

        Parameters
        ----------
        model : :class:`~compas_dem.models.SurfaceModel`
            The model to analyse
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_compl_energy_optimiser(
            printout=printout,
            plot=plot,
            max_iter=max_iter,
            starting_point=starting_point,
            support_displacement=support_displacement,
            Emethod=Emethod,
            solver=solver,
        )

        if printout:
            print("-" * 20)
            print("Complementary energy analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    @classmethod
    def create_lp_analysis(
        cls,
        formdiagram: FormDiagram,
        envelope: Envelope,
        solver: str = "CVXPY",
        printout: bool = False,
        plot: bool = False,
        max_iter: int = 500,
    ) -> "Analysis":
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`, optional
            The shape constraining the problem, by default None
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """
        optimiser = Optimiser.create_lp_optimiser(solver=solver, printout=printout, plot=plot, max_iter=max_iter)

        if printout:
            print("-" * 20)
            print("Load path Analysis created")
            print(optimiser)

        return cls(formdiagram, envelope, optimiser=optimiser)

    # =============================================================================
    # Methods
    # =============================================================================

    def is_convex(self) -> bool:
        """Check if the analysis problem is convex."""

        if not self.optimiser:
            raise ValueError("Define the Optimiser for the problem")

        objective = self.optimiser.settings["objective"]
        features = self.optimiser.settings["features"]

        if objective == "loadpath" and "fixed" in features:
            return True
        else:
            return False

    def set_optimiser_options(self, **kwargs):
        """Set the additional options of the optimisation."""
        if kwargs:
            self.optimiser.set_additional_options(**kwargs)

    def clear_previous_results(self):
        """Clear Previous results stored in the Analysis object.
        Necessary to perform sequential optimisation with the same analysis object"""
        self.optimiser.clear_optimiser()

    def apply_selfweight(self, normalize_loads=True):
        """Invoke method to apply selfweight to the nodes of the form diagram based on the envelope"""
        self.envelope.apply_selfweight_to_formdiagram(self.formdiagram, normalize=normalize_loads)

    def apply_hor_multiplier(self, multiplier=1.0, component="x"):
        """Apply a multiplier on the selfweight to the nodes of the form diagram based on the envelope

        Parameters
        ----------
        multiplier : float, optional
            Value of the horizontal multiplier, by default ``1.0``.
        component : str, optional
            Direction to apply the loads, ``x`` or ``y``, by default ``x``.

        Returns
        -------
        None
            The FormDiagram is modified in place.

        """
        arg = "p" + component

        for key in self.formdiagram.vertices():
            pz = self.formdiagram.vertex_attribute(key, "pz")
            self.formdiagram.vertex_attribute(key, arg, -1 * pz * multiplier)  # considers that swt (pz) is negative

    def apply_envelope(self):
        """Invoke method to apply ub and lb to the nodes based on the envelope's intrados and extrados"""

        self.envelope.apply_bounds_to_formdiagram(self.formdiagram)

    def apply_target(self):
        """Apply target to the nodes based on the shape's target surface"""

        self.envelope.apply_target_heights_to_formdiagram(self.formdiagram)

    def apply_envelope_on_xy(self, c=0.5):
        """_summary_

        Parameters
        ----------
        c : float, optional
            Distance in (x, y) to constraint the nodes limiting the hor. movement, by default 0.5

        """
        form = self.formdiagram
        for vertex in form.vertices():
            x, y = form.vertex_attributes(vertex, names=["x", "y"])

            form.vertex_attribute(vertex, name="xmin", value=x - c)
            form.vertex_attribute(vertex, name="xmax", value=x + c)
            form.vertex_attribute(vertex, name="ymin", value=y - c)
            form.vertex_attribute(vertex, name="ymax", value=y + c)

    def apply_cracks(self, key, position):
        """Apply cracks on the nodes (key) and in the positions up / down"""
        raise NotImplementedError

    def apply_reaction_bounds(self, assume_shape=None):
        """Apply limit thk to be respected by the anchor points"""
        self.envelope.apply_reaction_bounds_to_formdiagram(self.formdiagram)

    # =========================================================================
    # New Workflow Methods
    # =========================================================================

    def create_base_problem(self):
        """Create the base problem structure (lightweight, no solving).

        This method creates the fundamental problem structure from the form diagram
        and envelope, storing it in `optimiser.problem`. No optimization is performed.

        The base problem includes:
        - Basic equilibrium matrices (C, Ci, Cit)
        - Connectivity information
        - Load vectors
        - Bounds from envelope

        Notes
        -----
        This is the first step in the new workflow:
        1. create_base_problem()
        2. compute_starting_point()
        3. setup_optimization()
        4. solve()

        Examples
        --------
        >>> analysis = Analysis(form, envelope, optimiser)
        >>> analysis.create_base_problem()
        >>> print(f"Created problem with {len(analysis.optimiser.problem.edges)} edges")

        """
        set_up_base_problem(self)
        return self

    def compute_starting_point(self):
        """Compute and apply starting point to the form diagram.

        This method applies the starting point strategy specified in optimiser settings
        to initialize the form diagram. Strategies include:
        - "loadpath": Solve convex loadpath optimization
        - "tna": Use Thrust Network Analysis
        - "fdm": Use Force Density Method
        - "sag": Use sag approach
        - "current": Use current form state

        The computed starting point updates the form diagram's vertex z-coordinates
        and edge force densities (q).

        Notes
        -----
        This is the second step in the new workflow. It requires that
        `create_base_problem()` has been called first (for "loadpath" strategy).

        Examples
        --------
        >>> analysis.create_base_problem()
        >>> analysis.compute_starting_point()
        >>> print(f"Z range: {analysis.formdiagram.z_range()}")

        """
        # Local imports to avoid circular dependency
        from compas_tno.algorithms import equilibrium_fdm
        from compas_tno.solvers.startingpoint import startingpoint_loadpath
        from compas_tno.solvers.startingpoint import startingpoint_sag
        from compas_tno.solvers.startingpoint import startingpoint_tna

        form = self.formdiagram
        optimiser = self.optimiser
        problem = optimiser.problem
        settings = optimiser.settings
        starting_point = settings.get("starting_point", "current")
        boundary_force = settings.get("boundary_force", 50.0)
        printout = settings.get("printout_loadpath", True)
        find_inds = settings.get("find_inds", False)
        solver_convex = settings.get("solver_convex", "CLARABEL")

        if starting_point == "current":
            pass
        elif starting_point == "sag":
            startingpoint_sag(form, boundary_force=boundary_force)
        elif starting_point == "loadpath":
            startingpoint_loadpath(form, problem=problem, find_inds=find_inds, solver_convex=solver_convex, printout=printout)
        elif starting_point == "relax":
            equilibrium_fdm(form)
            startingpoint_tna(form)
        elif starting_point == "tna" or starting_point == "TNA":
            startingpoint_tna(form)
        else:
            print("Warning: define starting point")

        # Update problem.q from form after starting point is applied
        from numpy import array

        problem.q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)

        return self

    def setup_optimization(self):
        """Setup the full nonlinear optimization problem with all features.

        This method creates the complete optimization problem including:
        - All selected variables (q, zb, xyb, etc.)
        - All constraints (funicular, envelope, reactions, etc.)
        - Objective function and gradients
        - Jacobian matrices

        The problem is stored in `self.optimiser.problem` and includes all
        callable functions (fobj, fgrad, fconstr, fjac) ready for the solver.

        Notes
        -----
        This is the third step in the new workflow. It requires that:
        - Form diagram has been initialized with a starting point
        - Optimiser settings specify the features, variables, and constraints

        For convex problems, this method does nothing (base problem is sufficient).

        Examples
        --------
        >>> analysis.create_base_problem()
        >>> analysis.compute_starting_point()
        >>> analysis.setup_optimization()
        >>> print(f"Variables: {analysis.optimiser.problem.variables}")
        >>> print(f"Constraints: {analysis.optimiser.problem.constraints}")

        """
        if self.is_convex():
            # Convex problems don't need full optimization setup
            return self

        # Setup full nonlinear optimization problem
        set_up_general_optimisation(self)
        return self

    def solve(self, solver: Optional[str] = None, post_process: bool = True):
        """Solve the optimization problem and optionally post-process results.

        Parameters
        ----------
        solver : str, optional
            Solver to use ("IPOPT" or "SLSQP"). If None, uses the solver
            specified in optimiser settings. Default is None.
        post_process : bool, optional
            Whether to apply the solution to the form diagram after solving.
            Default is True.

        Returns
        -------
        SolverResult
            The optimization result containing xopt, fopt, success status, etc.

        Notes
        -----
        This is the fourth and final step in the new workflow. It requires that:
        - `setup_optimization()` has been called (for nonlinear)
        - OR `create_base_problem()` has been called (for convex)
        - Problem is fully configured in optimiser.problem

        The result is stored in `self.result`.

        Examples
        --------
        >>> analysis.create_base_problem()
        >>> analysis.compute_starting_point()
        >>> analysis.setup_optimization()
        >>> result = analysis.solve(solver="IPOPT")
        >>> print(f"Optimal objective: {result.fopt}")

        """
        # Determine solver
        if solver is None:
            solver = self.optimiser.settings.get("solver", "SLSQP")

        if not isinstance(solver, str):
            raise ValueError("Please provide the name of the solver")

        # Solve based on problem type
        if self.is_convex():
            result = run_convex_optimisation(self)
            self.result = result
        else:
            # Nonlinear optimization
            if solver.upper() == "IPOPT":
                result = run_nlopt_ipopt(self)
            elif solver.upper() == "SLSQP":
                result = run_nlopt_scipy(self)
            else:
                raise ValueError(f"Solver {solver} not supported. Please provide the an option between IPOPT or SLSQP")

            # Store the SolverResult in analysis
            self.result = result

        # Post-process if requested
        if post_process:
            post_process_nlopt(self)

        return self.result if hasattr(self, "result") and self.result else self

    # =========================================================================
    # Legacy Workflow Method
    # =========================================================================

    def set_up_optimiser(self):
        """With the data from the elements of the problem compute the matrices for the optimisation and the starting point"""
        # if self.is_convex():
        #     set_up_convex_optimisation(self)
        # else:
        #     set_up_general_optimisation(self)

        self.create_base_problem()
        self.compute_starting_point()
        self.setup_optimization()

    def run(self):
        """Run the complete optimization workflow (legacy convenience method).

        This method executes the full optimization workflow in one call:
        1. Sets up the problem (convex or nonlinear)
        2. Runs the optimization
        3. Post-processes the results

        Notes
        -----
        This is a convenience method that maintains backward compatibility.
        For more control over the workflow, use the new methods:
        - create_base_problem()
        - compute_starting_point()
        - setup_optimization()
        - solve()

        In a future version, this method will be refactored to call those
        methods internally instead of duplicating logic.

        """
        solver = self.optimiser.settings.get("solver", "SLSQP")

        if not isinstance(solver, str):
            raise ValueError("Please provide the name of the solver")

        if self.is_convex():
            result = run_convex_optimisation(self)
            self.result = result
            post_process_nlopt(self)
        else:
            if solver.upper() == "IPOPT":
                result = run_nlopt_ipopt(self)
            elif solver.upper() == "SLSQP":
                result = run_nlopt_scipy(self)
            else:
                raise ValueError(f"Solver {solver} not supported. Please provide the an option between IPOPT or SLSQP")

            # Store the SolverResult in analysis
            self.result = result

            post_process_nlopt(self)
