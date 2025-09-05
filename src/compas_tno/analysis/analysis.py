from typing import Annotated
from typing import Literal
from typing import Optional

import numpy.typing as npt

from compas.data import Data
from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import Envelope
from compas_tno.optimisers import Optimiser
from compas_tno.problems import set_up_convex_optimisation
from compas_tno.problems import set_up_general_optimisation
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

    def __init__(
        self,
        formdiagram: FormDiagram,
        envelope: Envelope,
        optimiser: Optimiser,
        settings: Optional[dict] = None,
        name: Optional[str] = None,
    ):
        name = name or "Analysis"

        super().__init__(name=name)

        self.settings = settings or {}
        self.formdiagram = formdiagram
        self.envelope = envelope
        self.optimiser = optimiser

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

    def set_up_optimiser(self):
        """With the data from the elements of the problem compute the matrices for the optimisation and the starting point"""
        if self.is_convex():
            set_up_convex_optimisation(self)
        else:
            set_up_general_optimisation(self)

    def run(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""
        solver = self.optimiser.settings.get("solver", "SLSQP")

        if not isinstance(solver, str):
            raise ValueError("Please provide the name of the solver")

        if self.is_convex():
            run_convex_optimisation(self)
        else:
            if solver.upper() == "IPOPT":
                self = run_nlopt_ipopt(self)
            elif solver.upper() == "SLSQP":
                self = run_nlopt_scipy(self)
            else:
                raise ValueError(f"Solver {solver} not supported. Please provide the an option between IPOPT or SLSQP")

            post_process_nlopt(self)
