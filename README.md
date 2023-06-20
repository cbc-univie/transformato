Transformato
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/wiederm/transformato/workflows/CI/badge.svg)](https://github.com/wiederm/transformato/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/wiederm/transformato/branch/master)
[![Github release](https://badgen.net/github/release/wiederm/transformato)](https://github.com/wiederm/transformato/)
[![GitHub license](https://img.shields.io/github/license/wiederm/transformato?color=green)](https://github.com/wiederm/transformato/blob/main/LICENSE)
[![GH Pages](https://github.com/wiederm/transformato/actions/workflows/build_page.yaml/badge.svg)](https://github.com/wiederm/transformato/actions/workflows/build_page.yaml)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/transformato.svg)](https://anaconda.org/conda-forge/transformato)

Transformato is a package that helps to set up relative alchemical free energy calculations of small molecules with a common core scaffold, either for solvation free energy[^1] or binding free energy[^2] estimates. The package is designed to be used with output generated by [CHARMM-GUI](https://charmm-gui.org/).

## Installation

transformato can be easily installed from conda-forge:
```
conda install transformato -c conda-forge
```

## Usage

For more information on how to use transformato please visit the [documentation](https://wiederm.github.io/transformato/)

## Maintainers

- Marcus Wieder <marcus.wieder@univie.ac.at> (University of Vienna)
- Johannes Karwounopoulos <johannes.karwounopoulos@univie.ac.at> (University of Vienna)


### Copyright

:copyright: 2022, Marcus Wieder, Johannes Karwounopoulos, Stefan Boresch


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[^1]:  Wieder, M., Fleck, M., Braunsfeld, B., Boresch, S., *J. Comput. Chem.* 2022, 1. [DOI](https://doi.org/10.1002/jcc.26877)
[^2]:  Karwounopoulos, J., Wieder, M., Boresch, S., *Front. Mol. Biosci.* 2022, 9, 954638 [DOI](https://doi.org/10.3389/fmolb.2022.954638
)
