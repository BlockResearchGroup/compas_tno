# COMPAS TNO

![build](https://github.com/BlockResearchGroup/compas_tno/workflows/build/badge.svg)
[![GitHub - License](https://img.shields.io/github/license/compas-dev/compas.svg)](https://github.com/BlockResearchGroup/compas_tno/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/COMPAS.svg)](https://pypi.python.org/project/COMPAS)
[![DOI](https://zenodo.org/badge/197799108.svg)](https://zenodo.org/badge/latestdoi/197799108)

> **⚠️ IMPORTANT NOTE** ⚠️
> 
> As of Summer 2025 COMPAS TNO is being refactored to work with COMPAS 2.0 and better integrate with the COMPAS MASONRY environment. A Legacy version is available at secondary branches.
> 
> **⚠️ IMPORTANT NOTE** ⚠️

**COMPAS TNO** is a Python package to find admissible thrust networks in masonry vaulted structures built in the [COMPAS](https://compas.dev/) framework.

Based on [Ricardo Maia Avelino](https://ricardoavelino.github.io/)'s doctoral thesis at ETH Zurich, this Package enables finding multi-objective particular internal stress solutions in masonry vaults, as the ones presented below.

![COMPAS TNO Objectives](./docs/_images/objectives.png)

## Installation

Stable releases can be installed from PyPI.

```bash
pip install compas_tno
```

To install the latest version for development, do:

```bash
git clone https://github.com/blockresearchgroup/compas_tno.git
cd compas_tno
pip install -e ".[dev]"
```

To install a version with support for IPOPT and Mosek

```bash
cd compas_tno
conda env create -f environment.yml
conda activate tno-dev
```

To get started with `compas_tno` have a look at [the documentation](https://github.com/blockresearchgroup/compas_tno).

## Questions and Feedback

For questions and feedback, have a look at the [COMPAS Forum](https://forum.compas-framework.org).

## Issues

If you run into problems, please file an issue on the [issue tracker](https://github.com/blockresearchgroup/compas_tno/issues). If we don't know it is broken, we can't fix it...

## Contributing

Guidelines for developers are under construction. However, we always accept contributions in the form of Pull Requests.

## Citing

If you use `compas_tno` for your research, cite one of our [papers](https://blockresearchgroup.github.io/compas_tno/latest/publications.html).

## License

`compas_tno` is licensed under the MIT License. See [LICENSE](https://github.com/blockresearchgroup/compas_tno/blob/main/LICENSE), for more information.
