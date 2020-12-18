
Cream User Guide
================

------------
Installation
------------

To install cream it, is best to use ``pip``. In folder ``cream``, ``setup.py``
file is provided that defines how to install cream using Python ``setuptools``.
Because, SCONE is intended as a prototyping environment it is best to set-up
cream package in so-called `editable` mode. This will allow ``pip`` to keep
track to any changes you make and keep the package up-to-date with the code
in the repository. To install go to the ``scone/cream`` directory and type::

    pip install -e .[test]

This will download all required dependencies. ``[test]`` is included to download
``pytest`` package, which is required to run unit tests.

Now you can verify that set-up was successful. First you can see if ``cream``
executable is in your environment path by typing::

    cream --help

You can also list all your packages handled by ``pip`` and find cream with::

    pip list | grep 'cream'

To keep your environment tidy, remember that you can remove cream by using::

    pip uninstall cream

This will just remove cream package from Python environment. Files in the
repository will not be removed.

Dependencies
------------
Must be pre-installed:
  * Python >=3.6
  * ``pip``

Will be installed by ``pip``:
  * ``click`` Python package for command line interfaces
  * ``pytest`` if ``[test]`` is specified

-----------------
Basic Information
-----------------

Cream command line interface is organised into a number of blocks.
Then each block contains number of subcommands that can perform different
operations. Help info that will be displayed is dependent on a current block.
For example if you type::

  cream --help

It will display all available blocks under `commands` e.g.::

  Usage: cream [OPTIONS] COMMAND [ARGS]...

    Command line utility for SCONE Monte Carlo Code.

  Options:
      --help  Show this message and exit.

  Commands:
      data  Manage nuclear data.

However if you enter a specific block information will change. For
``cream data --help``::

  Usage: cream data [OPTIONS] COMMAND [ARGS]...

    Manage nuclear data.

  Options:
    --help  Show this message and exit.

  Commands:
    relaxation  Generate atomic relaxation data library from ENDF.

Now ``Commands:`` lists available operations.

----------
Data Block
----------

Data block contains subcommands intended to help with the management of
nuclear data.

Relaxation
----------
Relaxation command allows to convert ENDF-formatted atomic relaxation data into
format that can be used by Serpent or SCONE. User must provide names of two
files that will contain transition (relaxation) data and electronic ground state
configurations respectively.

For example, if ENDF data is in `relax` folder with `.endf` extension,
the following command will perform conversion:: 

  cream data relaxation -o1 trans -o2 ground --format SCONE ./relax/*.endf
