.. SCONE documentation master file, created by
   sphinx-quickstart on Wed Mar  6 11:51:01 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is SCONE?
==============

SCONE stands for **S**\ tochastic **C**\ alculator **O**\ f **N**\ eutron
Transport **E**\ quation and it is a Monte Carlo particle transport code.
Unlike its more famous cousins (e.g. MCNP, Serpent, OpenMC) it does not strive
for the fastest execution and application in practical, production-level
calculations. Instead it aims to be a flexible code for academic use by
graduate students, who can use it to learn more about the detail of Monte Carlo
calculations and prototype new algorithms as part of their projects.

With that in mind SCONE aims to:

* Minimise time between an idea and prototype implementation on semi-realistic
  test case.
* Be easy to pick up by an engineering student.
* Allow enough flexibility to easily implement non-standard calculations.
* Offer reasonable performance, but prioritise flexibility over execution speed.

If your aim is to simply obtain standard results (criticality, power
distribution) for a reactor concept, SCONE might be a wrong tool for that
purpose and we would recommend looking at excellent software like Serpent or
OpenMC.

On the other hand, if you have an idea you want to prototype and see if it is
worth the considerable effort to implement it into the more optimised codes,
we hope that you will find our framework useful.

UK Note
-------

We would like to stress that we do not have a preferred pronunciation. Both the
*tone* and *gone* is perfectly fine. However, it is a **strong** opinion of this
author that the jam goes on top the cream.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Dictionary Input
   Unit Testing
   Style Guide



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
