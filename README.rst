*****
SCONE
*****

SCONE (**S**\ tochastic **C**\ alculator **O**\ f **N**\ eutron Transport **E**\ quation)
is an attempt to create an Object-Oriented framework for Monte Carlo particle transport
calculations. Its primary goal is to support neutron transport problems in nuclear reactor
physics. Object-oriented structure was selected to allow for better modularity and
(hopefully) ability to modify solution algorithms with virtually no knowledge about
the implementation details of large chunks of functionality (tallies, geometry, nuclear data
handling). Fortran 2008 was chosen over other languages (mostly C++) for its performance combined
with relative simplicity.


Getting Started
===============

Prerequisites
-------------
* Cmake (>=3.0)
* Fortran compiler, gfortran (>=6.3).
* LAPACK and BLAS Libraries
* pFUnit unit testing framework and Python interpreter
* UNIX-like environment (e.g. Linux, MacOS or Cygwin in Windows)

Documentation
-------------

Sphinx documentation is available in the ``doc`` folder. It is readable with
any reStructuredText (RST) viewer, but it is best to compile to html. To do that
install Sphinx with following command ::

    pip install --user -U sphinx

Then navigate to ``doc`` folder and compile ::

    cd ./doc
    make html

Documentation website will now be in ``./doc/_build/html``

Installation
------------

Instructions are available in Sphinx documentation
