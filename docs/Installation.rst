.. _installation:

Installation
============

Requirements
''''''''''''

Following is required to compile SCONE:

* Cmake (>=3.0)
* Fortran compiler, gfortran (>=6.3).
* LAPACK and BLAS Libraries
* pFUnit unit testing framework and Python interpreter
* UNIX-like environment (e.g. Linux, MacOS or Cygwin in Windows)

Getting gfortran
''''''''''''''''
To verify that you have gfortran available by typing::

    gfortran --version

If it isn't or its version is too old you will need to get it. If you have root
access to your machine you can your package manager to install gfortran. On
Debian/Ubuntu Linux a command like that will suffice::

   sudo apt-get install gfortran

On other operating systems it might be significantly different. You will need to
find information on how to use package manager in your Linux distribution.
Pre-compiled gfortran packages can be found
`here <https://gcc.gnu.org/wiki/GFortranBinaries>`_

Without administrator privileges you may want to compile GCC from source.
Of course that requires having some a C compiler.

#. Download the source code of GCC from one of the
   `mirrors <https://gcc.gnu.org/mirrors.html>`_

#. Extract the archive and use a provided script to download all prerequisites::

      tar -xf gcc-9.1.0.tar.gz
      cd gcc-9.1.0
      ./contrib/download_prerequisites

#. Now configure the compilation::

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
If you have root access to your machine use package manager to obtain latest
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

     git clone https://Your-Bitbucket-Name@bitbucket.org/Your-Bitbucket-Name/scone.git

#. Create build folder in the project directory(e.g. Debug)::

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

Run Tests to Verify
'''''''''''''''''''

If you compiled SCONE with tests enabled (you should by the way) you can now
verify that it works correctly by running the automated test suites. You
**must** to execute the following commands from ``scone`` directory. Some
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
SCONE with it. We are looking forward to F18, but it does not yet support the
object-oriented features of Fortran.

Obtaining Nuclear Data
''''''''''''''''''''''

You need to contact one of the developers and we will guide you personally. Soon
a Python script to generate appropriate library file from ACE files will be
added to repository.
