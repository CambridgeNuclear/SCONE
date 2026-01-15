.. _installation:

Installation
============

Requirements
''''''''''''

.. admonition:: Required

   Fortran Compiler
     Currently SCONE requires gfortran (>=8.3). Support for other compilers is pending.

   CMake
     CMake cross-platform build system is used to run all configuration scripts. Version (>=3.10)
     is required.


   LAPACK and BLAS Libraries
     For the moment, SCONE requires LAPACK and BLAS linear algebra libraries. However, they are
     not used extensively in the code and this dependency will become optional or will be removed.

   GNU/Linux operating system
     In principle any UNIX-like operating system should support SCONE. However people have
     experienced some significant run-time bugs when running on MacOS. Since we do not have
     an access to a Mackintosh we were not able to identify the problem yet.

.. admonition:: Optional

   pFUnit 4 test framework and Python interpreter
     Both the unit and the integration tests in SCONE use pFUnit framework. To run it requires a
     python interpreter. NOTE that version 4 (contrarily to the older 3.0) requires the use of
     gfortran version 8.3 or newer. 

Getting gfortran
''''''''''''''''
To verify that you have gfortran available by typing::

    gfortran --version

If you do not or its version is too old you will need to get it. If you have root
access to your machine you can your package manager to install gfortran. On
Debian/Ubuntu Linux a command like that will suffice::

   sudo apt-get install gfortran

On other operating systems it might different. You will need to
find information on how to use package manager in your Linux distribution.
Pre-compiled gfortran packages can be found
`here <https://gcc.gnu.org/wiki/GFortranBinaries>`_

Without administrator privileges you may want to compile GCC from source.
Of course that requires having a C compiler.

#. Download the source code of GCC from one of the
   `mirrors <https://gcc.gnu.org/mirrors.html>`_

#. Extract the archive and use a provided script to download all prerequisites::

      tar -xf gcc-9.1.0.tar.gz
      cd gcc-9.1.0
      ./contrib/download_prerequisites

#. Now configure the compilation. The command below is crude. We only set a `prefix` where
   gcc will be installed after successful compilation and select languages (frontends) we want to
   include. Documentation of the configure script is
   `available <https://gcc.gnu.org/install/configure.html>`_ ::

      ./configure --prefix=/path/to/install --enable-languages=c,c++,fortran

#. Providing that the configuration was successful, you can now start compiling
   GCC. It is a large code base and the process can take as much as few hours.
   To speed it up you can use parallel compilation with ``make -j 8`` assuming
   that you have 8 processors available. You can use ``nproc`` in console to
   check how many are available. When ready type::

      make -j8
      make install

#. Once your finished now you need to modify some of your environmental
   variables to allow OS to find the executables and relevant libraries. In your
   `.bashrc`` add the following lines (depending on your install directory)::

      export export PATH=/path/to/install/bin:$PATH
      export LD_LIBRARY_PATH=/path/to/install/lib64:$LD_LIBRARY_PATH

Getting CMake
'''''''''''''
If you have root access to your machine use package manager to obtain the latest
version of CMake. If you don't you can follow the instructions.

#. Download the installation shell script for Linux from the
   `website <https://cmake.org/download>`_ e.g. `cmake-3.17.0-rc1-Linux-x86_64.sh`.

#. Then you can install CMake by typing and following the instructions::

      bash ./cmake-3.17.0-rc1-Linux-x86_64.sh

#. Add the CMake to you ``PATH`` in ``.bashrc``::

      export PATH=/cmake/install/folder/bin:$PATH


Installing pFUnit
'''''''''''''''''
This is only required if the unit tests are to be build.

#. Make sure python can be invoked by a command ``python`` by typing::

     python --version

#. Download the pFUnit repository and enter the source code folder::

     git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
     cd pFUnit

#. Create a build folder and compile the code::

     mkdir build
     cd build
     cmake ./..
     make tests
     make install

#. Export environmental variables required by pFUnit::

     export F90=gfortran
     export F90_VENDOR=GNU


LAPACK and BLAS
'''''''''''''''
If you have root access it is best to install these with your package manager.
Follow the instructions only if you want to compile LAPACK and BLAS from source

#. Download a version of LAPACK from `official website
   <http://www.netlib.org/lapack/>`_.

#. In some directory on your filesystem extract the archive.

#. Configure compilation with cmake by typing::

     mkdir Build
     cd Build
     cmake ./..

#. If you don't have a root access on your machine or you want to install LAPACK
   to  a custom directory, use ccmake to change CMAKE_INSTALL_PREFIX. In Build
   directory type::

     ccmake ./..
     <Navigate to CMAKE_INSTALL_PREFIX and change it to your folder>
     Press [c] to configure
     Press [g] to generate and exit

#. Now compile LAPACK and install by typing::

     make
     make install


Compiling SCONE
'''''''''''''''

#. If you want to install with tests set PFUNIT_INSTALL environmental variable
   to directory in which pFUnit was installed. It may be worth adding the line
   to your ``.bashrc`` ::

     export PFUNIT_DIR=~/pFUnit/build/

#. If your LAPACK installation is not in default system directories use
   LAPACK_INSTALL enviromental variable to help CMAKE find the library. e.g. ::

     export LAPACK_INSTALL=~/LAPACK

#. Download the repository. Run the following commands::

     git clone https://github.com/CambridgeNuclear/SCONE

#. Create build folder in the project directory (e.g. Build)::

     cd ./scone
     mkdir Build

#. Generate makefile with CMake and compile::

     cmake -E chdir ./Build cmake ./..
     make -C Build

#. To switch off compilation of tests use the following commands::

     cmake -E chdir ./Build cmake ./.. -DBUILD_TESTS=OFF
     make -C Build

#. Note that you can use ccmake utility to modify avalible options and
   regenerate your make file just type the following into your terminal and
   follow the instructions::

     ccmake ./Build

.. admonition:: CMake options

   LTO
     Enable link-time optimisation. It allows the compiler to perform extra optimisations between
     different compilation units (modules in Fortran). It is crucial for performance in SCONE, since
     it enables inlining of small type-bound procedures. Set to `ON` by default. To disable::

       cmake .. -DLTO=OFF

   COVERAGE
     Collect code coverage information. If `ON` it allows to use `lcov` and `genhtml` to create
     an HTML coverage report. It is `OFF` by default. Enable with::

       cmake -DCOVERAGE=ON

   BUILD_TESTS
     Build unit and integration tests. It is `ON` by default. If enabled, the pFUnit must be
     installed and PFUNIT_INSTALL set. To disable tests::

       cmake -DBUILD_TESTS=OFF

   DEBUG
     Enable extra run-time checks available in the compiler. It is `OFF` by default. To enable::

       cmake -DDEBUG=ON


Run Tests to Verify
'''''''''''''''''''

If you compiled SCONE with tests enabled (you should by the way) you can now
verify that it works correctly by running the automated test suites. You
**must** execute the following commands from ``scone`` directory. Some
integration tests use files in ``IntegrationTestFiles`` and have hard-coded
relative paths. **Integration tests may fail if they are run from other
directory**. Run::

    ./Build/unitTests
    ./Build/integrationTests

This assume that ``Build`` is the build directory. If the tests were successful
that is great. If some of them failed it is troubling. Please open an Issue in
the online repository so we can try to resolve what is the problem. Provide at
least the following information:

#. Compiler Used (with version)
#. Operating System

Unfortunately we do not have access to Intel Fortran compiler so we cannot test
SCONE with it. We are planning to add support for Flang soon.

Obtaining Nuclear Data
''''''''''''''''''''''

SCONE requires ACE-formatted nuclear data. The JEFF-3.3 evaluation can be download from the
OACD NEA `website <https://www.oecd-nea.org/dbdata/jeff/jeff33/>`__. In addition SCONE requires
its own library file. An example of it is given in *IntegrationTestFiles/testLib*. Its format is::

  ! This is a comment line
  ! Each line needs to contain three entries
  ! ZAID   Line Number   PATH
  92233.03c;  1;       <absolute_path>/9233JEF33.ace;
  1001.03c;   4069;    <absolute_path>/1001JEF33.ace;
  ...

`Line Number` is the line in the file at which a particular data card begins. Each line cannot
contain more then a single entry. Each component must be delimited by a semi-colon.

To generate the library file from the collection of raw ACE files one can use the
``scripts/make_ace_lib.sh`` bash script. It can be run with the following command:

.. code-block:: bash

  ./scripts/make_ace_lib.sh <output> <mode> <extension_suffix> <path_to_base_folder>

To get extra help run the script without any arguments. In the line above, <output> is the file that
will be created; <mode> can be either ``CE`` for continuous energy neutron data cards and ``SAB`` for
thermal scattering S(α,β) cards; <extension_suffix> is the final part of the extension that is common to
all the ace files of interest; <path_to_base_folder> is the path to the folder than contains the ace 
files. The script can search recursively inside folder structures too. An example call is:

.. code-block:: bash

  ./scripts/make_ace_lib.sh ./endfb8.aceXS CE nc ./Lib80x/

Sadly the script can search only for a single type of card in one pass. Thus to create a full
library with thermal data we need to do the following:

.. code-block:: bash

  ./scripts/make_ace_lib.sh ./tempCE  CE  ace ./path_to_CE_ace_files/
  ./scripts/make_ace_lib.sh ./tempSAB SAB t   ./path_to_SAB_ace_files/
  cat tempCE tempSAB > fullLib.xsfile
