.. _installation:

Installation
============

Requirements
''''''''''''

.. admonition:: Required

   Fortran Compiler
     Currently SCONE requires gfortran (>=8.3). Support for other compilers is pending.

   CMake
     CMake is a cross-platform build system used to run all configuration scripts. Version (>=3.10)
     is required.

   LAPACK and BLAS Libraries
     For the moment, SCONE requires the linear algebra libraries LAPACK and BLAS. However, they are
     not extensively used in the code and this dependency will become optional or be removed.

   UNIX operating system
     In principles, SCONE is supported on any UNIX-like operating system. This includes Linux
     distributions such as Ubuntu or Debian, and macOS.

.. admonition:: Optional

   pFUnit 4 test framework and Python interpreter
     Both the unit and the integration tests in SCONE use the pFUnit framework, which requires a 
     python interpreter to be run. NOTE: version 4 (contrarily to the older 3.0) requires the use of
     gfortran version 8.3 or newer.

Linux distributions
-------------------

Installing gfortran
###################

Check that gfortran is available by typing::

    gfortran --version

If you do not have gfortran, or if its version is too old you will need to get/update it. If you
have root access to your machine you may use your package manager to install/update gfortran. 
On Debian/Ubuntu distributions a command like the one below will suffice::

   sudo apt-get install gfortran

Other operating systems may require different commands, and you may need to find information on how
to use the package manager in your specific Linux distribution. Pre-compiled gfortran packages can be 
found `here <https://gcc.gnu.org/wiki/GFortranBinaries>`_

Without administrator privileges you may need to compile GCC from source, which requires the use of a
C compiler.

#. Download the GCC source code from one of the `mirrors <https://gcc.gnu.org/mirrors.html>`_

#. Extract the archive and use the provided script to download all prerequisites::

      tar -xf gcc-9.1.0.tar.gz
      cd gcc-9.1.0
      ./contrib/download_prerequisites

#. Now configure the compilation. The below command is crude. We only set a `prefix` where
   gcc will be installed after successful compilation and select languages (frontends) we want to
   include. Documentation of the `configure`` script is available
   `here <https://gcc.gnu.org/install/configure.html>`_ ::

      ./configure --prefix=/path/to/install --enable-languages=c,c++,fortran

#. Provided that the configuration was successful, you can now compile GCC.
   Note that it is a large code base and the process can take up to a few hours;
   to speed it up you can use parallel compilation with ``make -j 8``, assuming
   that you have 8 processor cores available. You can use ``nproc`` in the console 
   to check how many cores are available. When ready, type::

      make -j8
      make install

#. Once compilation is complete, you will need to modify some of your environmental
   variables to allow the operating system to find the executables and relevant 
   libraries. In your ``.bashrc`` file, add the following lines (depending on your 
   install directory)::

      export export PATH=/path/to/install/bin:$PATH
      export LD_LIBRARY_PATH=/path/to/install/lib64:$LD_LIBRARY_PATH

Installing CMake
################

If you have root access to your machine use your package manager to obtain the latest
version of CMake. Else, you can follow the instructions below.

#. Download the installation shell script from the
   `website <https://cmake.org/download>`_ e.g. `cmake-3.17.0-rc1-Linux-x86_64.sh` for Linux.

#. Install CMake by typing and following the instructions::

      bash ./cmake-3.17.0-rc1-Linux-x86_64.sh

#. Add CMake to your ``PATH`` environmental variable in ``.bashrc``::

      export PATH=/cmake/install/folder/bin:$PATH

Installing pFUnit
#################

Note: the following is only required if unit tests are to be built.

#. Check that python can be invoked by typing::

     python --version

#. Download the pFUnit repository and enter the source code folder::

     git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
     cd pFUnit

#. Create a build folder (e.g. build) and compile the source code using CMake::

     mkdir build
     cd build
     cmake ./..
     make tests
     make install

#. Export environmental variables required by pFUnit in your ``.bashrc`` file::

     export F90=gfortran
     export F90_VENDOR=GNU

LAPACK and BLAS
###############

If you have root access it is recommended to install these with your package manager.
Follow the instructions below only if you want to compile LAPACK and BLAS from source.

#. Download a version of LAPACK from `official website
   <http://www.netlib.org/lapack/>`_ and extract the archive in some directory of your
   filesystem.

#. Create a build directory (e.g. Build) and configure the compilation with CMake by 
   typing::

     mkdir Build
     cd Build
     cmake ./..

#. If you don't have root access on your machine or want to install LAPACK
   to a custom directory, use ccmake to change CMAKE_INSTALL_PREFIX. In the Build
   directory type::

     ccmake ./..
     <Navigate to CMAKE_INSTALL_PREFIX and change it to your folder>
     Press [c] to configure
     Press [g] to generate and exit

#. Now compile LAPACK and install it by typing::

     make
     make install

macOS
-----

Note: the installation tutorial for macOS assumes that you have root access to your
machine and makes use of the `Homebrew` package manager; however, you may use a
different package manager (e.g. `Anaconda`) if you are more familiar with it.

#. Check that your Mac is running on macOS >= 15.0. You may check the version of your
   operating system and update it if necessary by going into *System Settings* > *General* 
   > *Software Update*.

#. Install `Xcode` from the App Store. `Xcode` contains crucial headers which are read 
   and interpreted when compiling software containing C/C++ languages. Once installed, 
   launch `Xcode` so that it can complete its initialisation. A dialog will be presented 
   indicating which Simulator runtimes are built-in, and which Simulator runtimes you may 
   download. Choose `Continue` to finish setting up `Xcode`.

#. Open a new `Terminal` window. If `Terminal` is not docked, you may find it by opening 
   a new `Finder` window, then going to *Applications* > *Utilities*.

#. Install `Homebrew` by typing the following command in your `Terminal` window::

	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#. Once `Homebrew` is installed, type the following command in your `Terminal` window. 
   This will install the latest versions of all the packages required to correctly set 
   up and run SCONE::

     brew install gcc cmake python git openblas lapack libomp

#. Close your `Terminal` window. Open a new `Finder` window and navigate to your `Home` directory 
   (⌘ + ⇧ + h). Display hidden files (⌘ + ⇧ + .) and find the ``.zprofile`` file (this is the macOS 
   equivalent of the ``.bashrc`` file on Linux distributions). Open it and insert **any of the 
   following lines which are not already present** (note: this depends on whether you have a Mac 
   running on an Intel CPU or an ARM - Apple Silicon - chip):
   
   * Intel::
   
          # Setting PATH for Python 3.13. The original version is saved in .zprofile.pysave.
          PATH="/Library/Frameworks/Python.framework/Versions/3.13/bin:${PATH}"
          export PATH

          # Set shell environment for Homebrew.
          eval "$(/usr/local/bin/brew shellenv)"

          # Export pFUnit installation folder.
          export PFUNIT_DIR=~/pFUnit/build/

          # Export environmental variables required by pFUnit.
          export F90=gfortran
          export F90_VENDOR=GNU

          # Export OpenMP root path and flags.
          export OpenMP_ROOT=$(brew --prefix)/opt/libomp
          export LDFLAGS="-L/usr/local/opt/libomp/lib"
          export CPPFLAGS="-I/usr/local/opt/libomp/include"

   * ARM::
   
          # Setting PATH for Python 3.13. The original version is saved in .zprofile.pysave.
          PATH="/Library/Frameworks/Python.framework/Versions/3.13/bin:${PATH}"
          export PATH

          # Set shell environment for Homebrew.
          eval "$(/opt/homebrew/bin/brew shellenv)"

          # Export pFUnit installation folder.
          export PFUNIT_DIR=~/pFUnit/build/

          # Export environmental variables required by pFUnit.
          export F90=gfortran
          export F90_VENDOR=GNU

          # Export OpenMP root path and flags.
          export OpenMP_ROOT=$(brew --prefix)/opt/libomp
          export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
          export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

#. Save the changes you made in your ``.zprofile`` file and close it. 
   You may now hide hidden files (⌘ + ⇧ + .).

#. Open a new `Terminal` window. By default, it should open in your `Home` directory, 
   but if not navigate to it by entering::

     cd

#. Download the pFUnit repository from Git, enter the source code repository and 
   create a build directory (e.g. build) by typing the following commands::

	git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
	cd pFUnit
	mkdir build
	cd build

#. Before proceeding, **make sure that the default C compiler is Apple Clang by entering 
   the following command**::

     gcc --version
   
   If it is not, then you have an alias (symlink) pointing to another C compiler. In this 
   case, you have two options:
   
   * Remove the alias, which will default the C compiler back to Apple Clang for all future 
     compilations. To do so, open a new `Finder` window then open the ‘Go to Folder’ prompt 
     by pressing (⇧ + ⌘ + g) and entering /usr. Navigate to /local/bin, locate the `gcc` 
     alias and delete it. Once this is done, you may revert to your `Terminal` window and type::
   
	gcc --version
    
     to ensure that the default C compiler is Apple Clang. Now initialise CMake (you should 
     still be in the build folder on your `Terminal`) by typing::

	cmake ./..

   * Initialise CMake by specifying which C compiler to use. In your `Terminal` window enter 
     the following::

	cmake -D CMAKE_C_COMPILER=CLANG ./..

#. Compile tests and install by typing::

	make tests
	make install

Compiling SCONE
---------------

#. If you want to install SCONE with unit tests, set the PFUNIT_INSTALL environmental 
   variable to the directory in which pFUnit was installed. It may be worth adding the
   following line to your ``.bashrc`` file::

     export PFUNIT_DIR=~/pFUnit/build/

#. If your LAPACK installation is not in default system directories use
   LAPACK_INSTALL enviromental variable to help CMake find the library, e.g. ::

     export LAPACK_INSTALL=~/LAPACK

#. Download the SCONE repository using Git by typing::

     git clone https://github.com/CambridgeNuclear/SCONE

#. Create a build folder (e.g. Build) in the project directory::

     cd ./scone
     mkdir Build

#. Generate a make file with CMake and compile the source code::

     cmake -E chdir ./Build cmake ./..
     make -C Build

#. To switch off tests compilation use the following commands::

     cmake -E chdir ./Build cmake ./.. -DBUILD_TESTS=OFF
     make -C Build

#. Note that you can use the ccmake utility to modify available options and
   regenerate your make file by typing the following into your terminal and
   following the instructions::

     ccmake ./Build

.. admonition:: CMake options

   LTO (Link-time optimisation)
     Allows the compiler to perform extra optimisations between different compilation units 
     (modules in Fortran). It is crucial for performance in SCONE, since it allows inlining 
     of small type-bound procedures. `ON` by default. To disable it, compile with::

       cmake .. -DLTO=OFF

   COVERAGE
     Collects code coverage information. Allows the use of `lcov` and `genhtml` to create an
     HTML coverage report if `ON`. `OFF` by default. To enable it, compile with::

       cmake -DCOVERAGE=ON

   BUILD_TESTS
     Builds unit and integration tests. Requires pFUnit to be installed and the PFUNIT_INSTALL
     environmental variable to be set. `ON` by default. To disable it, compile with::

       cmake -DBUILD_TESTS=OFF

   DEBUG
     Enables extra run-time checks available in the compiler. `OFF` by default. To enable it,
     compile with::

       cmake -DDEBUG=ON

Running automated tests
-----------------------

If tests were enabled during the compilation of SCONE (recommended), you may now 
verify that it correctly works by running the automated test suites. Note that some
integration tests use files in the ``IntegrationTestFiles`` directory and have 
hard-coded relative paths. **As such, you must execute the following commands 
from the** ``scone`` **directory. Integration tests may fail if they are run from 
other directories**::

    ./Build/unitTests
    ./Build/integrationTests

This assumes that ``Build`` is the build directory. If any of the tests fail, 
please open an issue `here <https://github.com/CambridgeNuclear/SCONE/issues>`_ 
so we can investigate the problem. Provide at least the following information:

#. Compiler used and version
#. Operating system

Unfortunately, we do not have access to Intel Fortran compilers so we cannot test
SCONE on them. We are planning to add support for Flang soon.

Obtaining nuclear data
----------------------

SCONE requires ACE-formatted nuclear data to run actual simulations. The necessary data can be
downloaded from the OACD NEA `website <https://www.oecd-nea.org/dbdata/jeff/jeff33/>`__. Please
make sure to download both the `Neutron` (293K) and `Neutron TSL` files in `ACE` format, and
extract the archives in some directory of your choice. In addition, SCONE requires its own library 
file, whose format is given below::

  ! This is a comment line
  ! Each line needs to contain three entries
  ! ZAID   Line Number   PATH
  92233.03c;  1;       <absolute_path>/9233JEF33.ace;
  1001.03c;   4069;    <absolute_path>/1001JEF33.ace;
  ...

Here, ``Line Number`` is the line in the file at which a particular data card begins. Each line cannot
contain more than one entry, and each component must be delimited by a ';'. An example of such a file 
is given in *IntegrationTestFiles/testLib*.

To generate the library file from the collection of downloaded raw ACE files, one can use the
``scripts/make_ace_lib.sh`` bash script, which can be run using the following command (to get extra 
help, run the script without any arguments):

.. code-block:: bash

  ./scripts/make_ace_lib.sh /path/lib.xsfile CE ./path_to_ace_files/*.ace

The ``CE`` letters allow to switch between searching for continuous energy (CE) and thermal scattering S(α,β) 
(SAB) neutron data cards. Sadly, the script can search only for a single card type in one pass; thus, to create 
a full library with thermal scattering data included one needs to run the following:

.. code-block:: bash

  ./scripts/make_ace_lib.sh ./tempCE CE ./path_to_CE_ace_files/*.ace
  ./scripts/make_ace_lib.sh ./tempSAB SAB ./path_to_SAB_ace_files/*.ace
  cat tempCE tempSAB > fullLib.xsfile

Running your first simulation with SCONE
----------------------------------------

Once the ``fullLib.xsfile`` has been generated, we can run our first actual simulation. SCONE uses text files as
simulation inputs. Instances of such files are included in the ``InputFiles`` directory. For this example, we will
use the `JEZ` input file. Navigate to the ``InputFiles`` directory, open the `JEZ` file and locate the following
lines::
     
     nuclearData {

          handles {
               ceData { type aceNeutronDatabase; aceLibrary $SCONE_ACE;}
          }

          ...

You may replace the ``$SCONE_ACE`` environmental variable in the `JEZ` file with the absolute path to the ``fullLib.xsfile`` 
file or, better yet, export this variable in your ``.bashrc``/``.zprofile`` file (depending on your OS) by adding the 
following line::
     
     export SCONE_ACE=path_to_fullLib.xsfile

and saving the changes in either case. Note that if you export this variable in your ``.bashrc``/``.zprofile`` file, you
will need to close and re-open your `Terminal` window to apply the changes. Once this is done, from your `Terminal` window
navigate to the ``SCONE/Build`` directory and run the following command::

     ./scone.out ../InputFiles/JEZ

which will run SCONE using the `JEZ` input file. Congratulations on running your first SCONE simulation!