# SCONE
[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](LICENCE)

SCONE (**S**tochastic **C**alculator **O**f **N**eutron Transport **E**quation) is an object-oriented Monte Carlo
particle transport code for reactor physics. It is intended as an accessible environment for
graduate students to test and develop their ideas before contributing them to more established
codes suitable for design calculations.

## Prerequisites
Required
- Cmake (>=3.0)
- Fortran compiler, gfortran (>=6.3)
- LAPACK and BLAS Libraries
- GNU/Linux operating system

Optional
- pFUnit test framework
- Python 3 interpreter

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
