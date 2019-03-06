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

SCONE itself
''''''''''''
#. If you want to install with tests set PFUNIT_INSTALL environmental variable to directory in 
   which pFUnit was installed. e.g. :: 
   
     export PFUNIT_INSTALL=~/pFUnit    

#. If PFUNIT_INSTALL is unset and tests are on, cmake will attempt to download, compile and install 
   pFUnit into "external" folder in cmake build directory. This may take a long time so it is 
   recommended that pFUnit is installed beforehand and environmental variable is set.  

#. Download the repositry. Run the following commands:: 

     git clone https://Mikolaj_Adam_Kowalski@bitbucket.org/Mikolaj_Adam_Kowalski/cued-mc-code.git  
    
#. Create build folder in the project directory(e.g. Debug):: 

     cd ./cued-mc-code
     mkdir Debug
   
#. Generate makefile with CMake and compile::

     cmake -E chdir ./Debug cmake ./..
     make -C Debug

#. To switch off compilation of tests use the following commands:: 

     cmake -E chdir ./Debug cmake ./.. -DBUILD_TESTS=OFF 
     make -C Debug 

#. Note that to change form compilation with tests to without tests you need to delete all files
   in build directory (Debug in the example) and execute appropriate cmake command again. 
   
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
 
Style Guide (Draft)
===================
Please note that this style guide is new and a lot of code was written without it. As such it is 
likely that some portions of the code do not conform to following rules. If found, such section should 
be brought to my attention (mak60@cam.ac.uk). However, feel free to correct it yourself.     


Rules on files
--------------
#. Each source file must contain only a single module or program. All procedures (subroutines & 
   functions) must be contained within modules. No old fortran style separate files for procedures 
   are allowed.   

#. All source files should be free from and have a ``.f90`` extension

#. It is preferable to have a single public class per module. It is not a strict rule and sometimes 
   it is necessary to group number of classes together to avoid circular module dependencies.However, 
   such cases should be avoided and clearly marked in module description. Also name of the module should 
   correspond to name of the class. E.g. module ``myType_class`` must contain type ``myType``.

#. All files containing programs, should be placed into "Apps" folder of the source directory. No 
   programs outside "Apps" folder are allowed. 

#. Module file name must have exactly the same name as the module it contains.  
   
#. Content of a module file must be clearly identified by its suffix
     * ``_func`` contains a library of functions. (genericProcedures.f90 breaks it for now) 
     
     * ``_inter`` contains and abstract class (abstract interface). 
     
     * ``_class`` contains a class. 
     
     * Module with no suffix contains only global parameters (constants) 

Rules on code
-------------
#. Every source file needs to use ``implicit none``. No implicit typing is allowed. If you don't know 
   what implicit typing is, don't think about it. Just make shure that implicit none is present at 
   the beginning of a module or program e.g. :: 
   
     module properModule_func 
       implicit none 
       private 
       ! The above switches off implicit typing in the whole module 
       ! including module procedures. 
     
       public :: printHello
     
     contains 
       
       !!
       !! Prints greating from Fortran. 
       !!
       subroutine printHello() 
         ! No "implicit none" is needed here. 
          
         print *, "Hello, I am Fortran. I am a good programing language!"
          
       end subroutine printHello 
     
     end module properModule_func     

#. Keep variable, classes and procedure names descriptive. Try to keep them short. 
   In general use lowerCamelCase, but break this rule for short variables like 
   ``N`` for an integer or ``V`` for velocity or any other physical variable with obvious meaning. 
   
#. Always include ``numPrecision`` module in a source file and use parametrisation of the variables e.g. :: 
     
     program myProg 
       use numPrecision 
       implicit none 
       
       integer(shortInt)  :: i           ! standard integer 
       integer(longInt)   :: score       ! integer where large values are expected 
       real(defReal)      :: float       ! Every real number 
       logical(defBool)   :: isHappy     ! Every logical variable 
       character(nameLen) :: neutronName ! Use for names of of objects, short strings 
       character(pathLen) :: filePath    ! Use to store paths to file  
     
     end program myProg          

#. Always explicitly import variables, types and functions from modules, unless module has now suffix 
   and contains only global parameters. :: 
   
     program prog 
       use numPrecision 
       use endfConstants 
       use dictionary_class,  only : dictionary 
       use nuclearData_inter, only : nuclearData 
       ...
     end program prog
      
.. Finish stuff beyond this point 

#. Every procedure definition needs to specify intent for its dummy variables. Type of the function 
   is to be defined in its variable declarations. Dummy arguments must be defined in order of their 
   apperance in argument list. Result type should be declared immediately after dummy arguments. Local 
   variables are to be defined next.   :: 
     
     pure recursive function factorial(n) result(fact) 
       integer(shortInt), intent(in) :: n    ! Value n connot be modified (attempt will produce compiler error) 
       integer(shortInt)             :: fact ! Define type of function result 
        
        
        
     end function factorial 

#. Every type or procedure needs to contain comment above itself marked with ``!!`` and with short 
   description of the prodedure or the type.
     
