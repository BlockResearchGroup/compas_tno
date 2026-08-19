# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added specialized Problem classes: `FixedProblem`, `SymmetricProblem`, `FixedSymmetricProblem` with clean `from_formdiagram()` constructors.
* Added class-based solver API: `IPOPTSolver`, `ScipySolver`, `CVXPYSolver` inheriting from abstract base classes.
* Added `SolverResult` dataclass with standardized attributes: `xopt`, `fopt`, `success`, `message`, `niter`, `time`, `exitflag`.
* Added `apply_solution_to_form()` function for decoupled post-processing using `SolverResult`.
* Added helper method `Problem._to_specialized()` to eliminate code duplication when creating specialized problem types.
* Added helper functions in `setup.py` to improve readability: `_apply_bounds_to_edges()`, `_create_problem_from_features()`, `_setup_problem_metadata()`, `_apply_starting_point()`, `_build_variable_vector()`, `_compute_initial_values()`, `_print_diagnostics()`.

### Changed

* Refactored `Problem` class to store optimization functions directly: `fobj`, `fconstr`, `fgrad`, `fjac` (callable functions) distinct from `variables` and `constraints` (lists of strings).
* Refactored `Problem` class to store optimization variables: `x0`, `bounds`, `f0`, `g0` for simplified solver interface.
* Refactored `set_up_general_optimisation()` from 363 lines to ~116 lines using extracted helper functions.
* Refactored `FixedSymmetricProblem` to build on `FixedProblem` instead of duplicating logic, eliminating 41 lines of code duplication.
* Updated solver classes to automatically retrieve functions from `problem` if not explicitly provided, enabling simplified API: `solver.solve(problem)`.
* Updated `NonlinearSolver` base class to use stored `problem.x0` and `problem.bounds` when available.
* Updated `post_process_nlopt()` to use `apply_solution_to_form()` internally with `SolverResult` object.
* Updated tests to use new clean API exclusively.
* **Architectural improvement**: Results now stored in `analysis.result` as `SolverResult` object instead of scattered across `optimiser` attributes.
* Moved `startingpoint` module from `problems/` to `solvers/` to avoid circular imports and improve logical organization.

### Removed

* Removed deprecated constructor functions: `initialise_form()`, `initialise_problem_general()`.
* Removed deprecated adaptor functions: `adapt_problem_to_fixed_diagram()`, `adapt_problem_to_sym_diagram()`, `adapt_problem_to_sym_and_fixed_diagram()`.
* Removed 430 lines of duplicated code from `problems.py` (33% reduction) through refactoring.


## [0.3.0] 2025-09-06

### Added

* Separated routines for `starting_point`.

### Changed

* Updated `Analysis` object and optimization routines to work with the new infrastructure.

### Removed

* Removed modules `.diagrams`, `.shape`. TNO is meant to operate with `compas_tna` classes. 
* Removed module`.viewer`. New viewer is in `compas_masonnry`.
* Removed `MATLAB` convex optimization to simplify installation
* Removed unused solvers such as `MMA` and `PyOpt`. New focus is on `slsqp` and `ipopt`.
* Made `ipopt` optional to simplify base installation.

## [0.2.2] 2023-09-03

### Added

- Added new functions in FormArtists
- Corrected typo in function names
- Improved sync with Proxy
- Added QR decomposition in solvers

### Changed

### Removed


## [0.2.1] 2023-05-04

### Added

- Updated workflow.

### Changed

### Removed

## [0.2.0] - 2023-05-03

### Added

- Initial version.

### Changed

### Removed
