
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

#. All source files should be free from and have a ``.F90`` extension. Note that the capital letter is 
   a convention that tells the compiler to use a preprocessor on the given file.  

#. It is preferable to have a single public class per module. It is not a strict rule and sometimes 
   it is necessary to group number of classes together to avoid circular module dependencies. However, 
   such cases should be avoided and clearly marked in module description. Also name of the module should 
   correspond to name of the class. E.g. module ``myType_class`` must contain type ``myType``.

#. All files containing programs, should be placed into "Apps" folder of the source directory. No 
   programs outside "Apps" folder are allowed. 

#. Module file name must have exactly the same name as the module it contains.  
   
#. Content of a module file must be clearly identified by its suffix:
     * ``_func`` contains a library of functions. (genericProcedures.f90 breaks it for now) 
     
     * ``_inter`` contains and abstract class (abstract interface). 
     
     * ``_class`` contains a class. 
     
     * ``_mod`` indicates a class-like (object-like) module. It is effectively an implementation of 
       a Singleton design pattern. The module contains some data and access function to manipulate it. 
       Its state can change during exectution.   
     
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

#. We follow the rule that source code should not extend beyond the column 100. 
   
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

#. If you need to import a lot of modules group the use statements in some logicaly orderd blocks with  
   some associated comment. Try to make reading of imported modules easy :: 
     
     program prog 
       use numPrecision 
       use endfConstants 
       use dictionary_class,   only : dictionary 
       use nuclearData_inter,  only : nuclearData 
       
       ! Import factories 
       use myInterface_func,   only : new_myInterface 
       use yourInterface_func, only : new_yourInterface
       ...             
     end program prog
     
#. Every procedure definition needs to specify intent for its dummy variables. Type of the function 
   is to be defined in its variable declarations. Dummy arguments must be defined in order of their 
   apperance in argument list. Result type should be declared immediately after dummy arguments. Local 
   variables are to be defined next. If procedure can return errors it should have a character parameter
   Here witch contains name of the function and file in which it is defined. :: 
     
     pure recursive function factorial(n) result(fact) 
       integer(shortInt), intent(in) :: n    ! Value n connot be modified (attempt will produce compiler error) 
       integer(shortInt)             :: fact ! Define type of function result 
       integer(shortInt)             :: temp ! Define local variable 
       character(100), parameter :: Here ='factorial (math_func.F90)' ! Location information for error reporting    
       ... 
     end function factorial 

#. Every type or procedure needs to contain comment above itself marked with ``!!`` and with short 
   description of the prodedure or the type. 
   
#. Try to follow the following pattern for the procedure description ::
   
     !!
     !! Brief description of what myFunc does 
     !!
     !! Detailed Description of what myFunc does 
     !!
     !! Args: 
     !!   A [intent] -> explenation of the argument 
     !!   B [intent] -> exlpenation of the argument may be long so 
     !!     it is necessary to move it the the next line sometimes 
     !!
     !! Result:
     !!   Describe what the result of the procedure is.    
     !!                
     !! Errors:
     !!   Describe how does the procedure behaves for invalid arguments as well as under what 
     !!   conditions it fails 
     !!
     !! NOTE: Any important information you want to highlight 
     !!  
     function myFunc(A,B) result(C) 
       Procedure Definition  
     end function 

#. Try to follow the following pattern for the derived type(clas) description. Use the same format 
   for the class-like modules (with _mod suffix) :: 
     
     !!
     !! Brief description of the type 
     !!
     !! Detailed description of the type 
     !!
     !! Public Members: 
     !!   A -> Brief description of memebers. May be missing if there are none 
     !!
     !! Private Memebers:
     !!   B -> Brief Description of members. May be missing if there are none 
     !!
     !! Interface: 
     !!   method1 -> Brief Description of class methods. Very short. Details should 
     !!     be contained within comment above procedure definition.  
     !!
     type myType
       Type Definition 
     end type    