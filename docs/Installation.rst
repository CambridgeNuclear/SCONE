.. _installation:

Installation
============

Requirements
''''''''''''

.. admonition:: Required

   Fortran Compiler
     Currently SCONE requires gfortran (>=6.3). Support for other compilers is pending.

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

   pFUnit 3 test framework and Python interpreter
     Both the unit and the integration tests in SCONE use pFUnit framework. To run it requires a
     python interpreter. Not that we use version 3.0 despite, newer 4.0 being available. This is
     to retain support for gfortran in versions older then 8.3 (required by pFUnit 4.0).

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

#. Create a folder for the local installation of pFUnit e.g. in your home
   directory and download the pFUnit repository and enter the source code folder::

     mkdir pFUnit
     cd pFUnit
     git clone git://git.code.sf.net/p/pfunit/code pfunit-code
     cd pfunit-code

#. Export environmental variables required by pFUnit::

     export F90=gfortran
     export F90_VENDOR=GNU

#. Build and test pFUnit by typing::

     make tests

#. Install pFUnit in any directory you have access to e.g. ::

     make install INSTALL_DIR=~/pFUnit

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

     export PFUNIT_INSTALL=~/pFUnit

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
  92233.03c  1       <absolute_path>/9233JEF33.ace
  1001.03c   4069    <absolute_path>/1001JEF33.ace
  ...

`Line Number` is the line in the file at which a particular data card begins. Each line cannot
contain more then a single entry. Each component must be space separated. Path can have only 100
character and contain no spaces.

Soon the format of the file will be changes to allow spaces in the path. Also the limitation on its
length will be lifted. A script will be included in ``cream`` to generate the library file directly
from the ACE files.
