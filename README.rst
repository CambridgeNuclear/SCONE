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

* Modern Fortran Compiler 

  * gfortran 6.3 or higher 
* CMake 3.0.0 or higher 
* LAPACK and BLAS Library 
* pFUnit framework if testing is to be enabled

  * pFUnit requires python to be installed  


Installation
------------

Unit Testing Framework
''''''''''''''''''''''
#. Make sure python can be invoked by a command 'python' by typing:: 

     python --version 

#. Create a folder for the local installation of pFUnit e.g. in your home directory and 
   download the pFUnit repository and enter the source code folder:: 
   
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

LAPACK AND BLAS
'''''''''''''''
#. Download a version of LAPACK from `official website <http://www.netlib.org/lapack/>`_.

#. In some directory on your filesystem extract the archive.

#. Configure compilation with cmake by typing:: 

     mkdir Build 
     cd Build
     cmake ./..

#. If you don't have a root access on your machine or you want to install LAPACK to  a custom 
   directory, use ccmake to change CMAKE_INSTALL_PREFIX. In Build directory type::
   
     ccmake ./..  
     <Navigate to CMAKE_INSTALL_PREFIX and change it to your folder> 
     Press [c] to configure 
     Press [g] to generate and exit 
     
#. Now compile LAPACK and install by typing:: 

     make 
     make install 
     
Compiling SCONE
'''''''''''''''

#. If you want to install with tests set PFUNIT_INSTALL environmental variable to directory in 
   which pFUnit was installed. e.g. :: 
   
     export PFUNIT_INSTALL=~/pFUnit    

#. If your LAPACK installation is not in default system directories use LAPACK_INSTALL enviromental 
   variable to help CMAKE find the library. e.g. :: 
   
     export LAPACK_INSTALL=~/LAPACK 

#. Download the repositry. Run the following commands:: 

     git clone https://Mikolaj_Adam_Kowalski@bitbucket.org/Mikolaj_Adam_Kowalski/scone.git  
    
#. Create build folder in the project directory(e.g. Debug):: 

     cd ./cued-mc-code
     mkdir Debug
   
#. Generate makefile with CMake and compile::

     cmake -E chdir ./Debug cmake ./..
     make -C Debug

#. To switch off compilation of tests use the following commands:: 

     cmake -E chdir ./Debug cmake ./.. -DBUILD_TESTS=OFF 
     make -C Debug 

#. Note that you can use ccmake utility to modify avalible options and regenerate your make file just 
type the following into your terminal and follow the instructions:: 

     ccmake ./Debug     

   
Development Guide (Draft)
=========================

Writing Tests in SCONE
----------------------
There are three main types of tests

* Unit Tests, which main goal is to enforce specification (interface and behaviour) of the code 
  component. In Fortran this component may be a procedure, derived type or a module. <Finish this 
  at some point> 

* Integration tests, which aim to verify that large pieces of the code behave as expected in
  **realistic use scenarios**. So for example particle crossing number of surfaces in PWR geometry 
  ends up in a place its supposed to.      
 
* Regression tests, which check that the output of the code for a given input has not changed.
  This should be a collection of input files for a real-life use cases. In Monte Carlo there is 
  a problem that if the sequence of the psudo-random number generator is changed regression tests
   will fail even if the result is still valid.      
 
<Finish by giving rules on how to write unit tests> 
 

     
