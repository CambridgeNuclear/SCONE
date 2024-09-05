# SCONE
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![build-and-test-ubuntu](https://github.com/CambridgeNuclear/SCONE/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/CambridgeNuclear/SCONE/actions/workflows/build-and-test.yml)
[![Documentation Status](https://readthedocs.org/projects/scone/badge/?version=latest)](https://scone.readthedocs.io/en/latest/?badge=latest)

SCONE (**S**tochastic **C**alculator **O**f **N**eutron Transport **E**quation) is an object-oriented Monte Carlo
particle transport code for reactor physics. It is intended as an accessible environment for
graduate students to test and develop their ideas before contributing them to more established
codes suitable for design calculations.

SCONE documentation is hosted at: <https://scone.readthedocs.io>

## Prerequisites
Required

* Cmake (>=3.10)
* Fortran compiler, gfortran (>=8.3)
* LAPACK and BLAS Libraries
* GNU/Linux operating system

Optional

* pFUnit 4 test framework
* Python 3 interpreter

## Installation
Instructions are avaliable in the Sphinx documentation.

## Compiling Documentation
Sphinx documentation is available in the docs folder. It is readable with any reStructuredText (RST)
viewer, but it is best to compile to html.

Compiling documentation requires few python packages. You can install them all with the following
command. Option `--user` installs them in your home directory and does not require administrator access.
```
pip install --user -U sphinx, sphinx_rtd_theme
```
Then natigate to `docs` folder and compile using `make`
```
make html
```

HTML documentation should now be avaliable in `./_build/html`

## Licence
This project is licensed under MIT Licence - see the [LICENCE](LICENCE) file for details.
