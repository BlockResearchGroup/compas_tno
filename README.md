# COMPAS TNO

![build](https://github.com/BlockResearchGroup/compas_tno/workflows/build/badge.svg)
[![GitHub - License](https://img.shields.io/github/license/compas-dev/compas.svg)](https://github.com/BlockResearchGroup/compas_tno/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/COMPAS.svg)](https://pypi.python.org/project/COMPAS)
[![DOI](https://zenodo.org/badge/197799108.svg)](https://zenodo.org/badge/latestdoi/197799108)

**COMPAS TNO** is a Python package to find admissible thrust networks in masonry vaulted structures built in the [COMPAS](https://compas.dev/) framework.

Based on [Ricardo Maia Avelino](https://ricardoavelino.github.io/)'s doctoral thesis at ETH Zurich, this Package enables finding multi-objective particular internal stress solutions in masonry vaults, as the ones presented below.

![COMPAS TNO Objectives](./docs/_images/objectives.png)

## Installation

The recommended an editable install of **COMPAS TNO** with [Anaconda/conda](https://conda.io/docs/). Here we create an environment called `tno` and install it:

```
conda create -n tno -c conda-forge python COMPAS triangle compas_view2
conda activate tno
git clone https://github.com/BlockResearchGroup/compas_tno.git
cd compas_tno
pip install -e .
```

## First Steps: Read the docs

A walkthrough the package is available in the documentation:
* <https://blockresearchgroup.github.io/compas_tno/>

## Issue tracker

If you find a bug, please help us solve it by [filing a report](https://github.com/BlockResearchGroup/compas_tno/issues).

## Citing

If you use **COMPAS TNO** for your research, cite one of our [papers](https://blockresearchgroup.github.io/compas_tno/latest/publications.html).

## License

**COMPAS TNO** is [released under the MIT license](https://github.com/BlockResearchGroup/compas_tno/latest/license.html).
