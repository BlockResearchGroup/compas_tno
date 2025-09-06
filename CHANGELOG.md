# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
