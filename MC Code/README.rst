*****
SCONE
*****

SCONE (**S**\ tochastic **C**\ alculations **O**\ f **N**\ eutron Transport **E**\ quation) 
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
* CMake 2.8.11.2 or higher 


Installation
------------

1. Download the repositry. Run the following commands:: 

     git clone https://Mikolaj_Adam_Kowalski@bitbucket.org/Mikolaj_Adam_Kowalski/cued-mc-code.git  
    
2. Create build folder in the project directory(e.g. Debug):: 

     cd ./cued-mc-code
     mkdir Debug
   
3. Generate makefile with CMake and compile::

     cmake -E chdir ./Debug cmake ./..
     make -C Debug



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
     
     * ``_inte`` contains and abstract class (abstract interface). 
     
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
     
