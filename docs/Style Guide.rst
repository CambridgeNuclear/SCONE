.. _style-guide:

Style Guide
===========
Please note that this style guide is new and a lot of code was written without
it. As such it is likely that some portions of the code do not conform to
following rules. If found, post their location in the `Doc Gaps` Issue, but feel
free to correct them yourself.


Rules on files
--------------
#. Each source file must contain only a single module or program. All procedures
   (subroutines & functions) must be contained within modules. No old FORTRAN
   style separate files for procedures are allowed.

#. All source files should be free from and have a ``.f90`` of ``.F90``
   extension. Note that the capital letter is a convention that tells the
   compiler to use a preprocessor on the given file.

#. It is preferable to have a single public class per module. It is not a strict
   rule and sometimes it is necessary to group number of classes together to
   avoid circular module dependencies. However, such cases should be avoided and
   clearly marked in module description. Also name of the module should
   correspond to name of the class. E.g. module ``myType_class`` must contain
   type ``myType``.

#. All files containing programs, should be placed into "Apps" folder of the
   source directory. No programs outside "Apps" folder are allowed.

#. Module file name must have exactly the same name as the module it contains.

#. Content of a module file must be clearly identified by its suffix

     * ``_func`` contains a library of functions. (genericProcedures.f90 breaks it for now)

     * ``_inter`` contains and abstract class (abstract interface).

     * ``_class`` contains a class.

     * ``_mod`` indicates a class-like (object-like) module. It is effectively
       an implementation of a Singleton design pattern. The module contains some
       data and access function to manipulate it. Its state can change during
       execution.

     * Module with no suffix contains only global parameters (constants)

Rules on code
-------------
#. Every source file needs to use ``implicit none``. No implicit typing is
   allowed. If you don't know what implicit typing is, don't think about it.
   Just make sure that implicit none is present at the beginning of a module or
   program e.g. ::

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
   ``N`` for an integer or ``V`` for velocity or any other physical variable
   with obvious meaning.

#. Always include ``numPrecision`` module in a source file and use
   parametrisation of the variables e.g. ::

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

#. Always explicitly import variables, types and functions from modules,
   unless module has now suffix and contains only global parameters. ::

     program prog
       use numPrecision
       use endfConstants
       use dictionary_class,  only : dictionary
       use nuclearData_inter, only : nuclearData
       ...
     end program prog

#. If you need to import a lot of modules group the use statements in some
   logically ordered blocks with some associated comment. Try to make reading
   of imported modules easy ::

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

#. Every procedure definition needs to specify intent for its dummy variables.
   Type of the function is to be defined in its variable declarations. Dummy
   arguments must be defined in order of their appearance in argument list.
   Result type should be declared immediately after dummy arguments. Local
   variables are to be defined next. If procedure can return errors it should
   have a character parameter ``Here``, which contains name of the function and
   file in which it is defined. ::

     pure recursive function factorial(n) result(fact)
       integer(shortInt), intent(in) :: n    ! Value n connot be modified (attempt will produce compiler error)
       integer(shortInt)             :: fact ! Define type of function result
       integer(shortInt)             :: temp ! Define local variable
       character(100), parameter :: Here ='factorial (math_func.F90)' ! Location information for error reporting
       ...
     end function factorial

#. When using numeric constants include correct parametrisation. ::

      real(defReal) :: x
      x = 9.76_defReal - OK
      x = 9.76_8       - WRONG
      x = 9.76         - WRONG

#. Note that ``numPrecision`` contains useful numeric parameters. ::

     real(defReal) :: x
     x = ONE
     x = TWO
     x = PI * TWO_PI

#. Use parameters in favour of numeric constants. They improve readability.
   Parameters should be in CAPITAL_LETTERS ::

    integer(shortInt), parameter :: X_AXIS = 1
    if(axis == X_AXIS) then    - OK
    if(axis == 1) then         - WRONG

Whitespaces and Indentation
---------------------------
#. Use no tabs. Only spaces. 2 Spaces per indentation level. It may be useful to
   configure your editor to insert 2 Spaces on Tab press.

#. We follow the rule that source code should not extend beyond the column 100.

#. When accessing components of derived types put space between % e.g. ::

     a = myType % type2 % func()  - OK
     a = myType%type2%func()      - WRONG

#. In procedure calls add a single space between arguments e.g. ::

     call mySubroutine(x, y(4), z)    - OK
     call mySubroutine( x, y( 2 ), z) - WRONG
     call mySubroutine(x,y(2),z)      - WRONG

#. It is useful to put whitespace around logical operators::

      if (x /= 2) then - OK
      if (x/=2) then   - WRONG


Comments and documentation
--------------------------

#. Every type or procedure needs to contain comment above itself marked with
   ``!!`` and with the description of the procedure or the type.

#. Follow the following pattern for the procedure description ::

     !!
     !! Brief description of what myFunc does
     !!
     !! Detailed Description of what myFunc does
     !!
     !! Args:
     !!   A [intent] -> explanation of the argument
     !!   B [intent] -> explanation of the argument may be long so
     !!     it is necessary to move it the the next line sometimes
     !!
     !! Result:
     !!   Describe what the result of the procedure is.
     !!
     !! Errors:
     !!   Describe how does the procedure behaves for invalid arguments as well
     !!   as under what conditions it fails. Describe only errors from execution
     !!   of this function. DO NOT include errors that may appear in procedures
     !!   called by the function.
     !!
     !! NOTE: Any important information you want to highlight
     !!
     function myFunc(A,B) result(C)
       Procedure Definition
     end function

#. Note that when giving errors information, errors from procedures called by
   the procedure we document should not in general be included in the *Errors:*
   section. This is because these sub-procedures may change and any changes to
   their error behaviour would (most likely) not be propagated to the
   documentation of all procedures that use them, rendering their description
   invalid. However it would be best to use common sense and indicate when
   a procedure is depending on errors given by the other procedure. Just try to
   make *Errors:* section informative. Use your best judgement. Usually this
   section is the most useful part of documentation so help you fellow users
   and your future self.

#. Try to follow the following pattern for the derived type(class) description.
   Use the same format for the class-like modules (with _mod suffix) ::

     !!
     !! Brief description of the type
     !!
     !! Detailed description of the type
     !!
     !! Public Members:
     !!   A -> Brief description of members. May be missing if there are none
     !!
     !! Private Members:
     !!   B -> Brief Description of members. May be missing if there are none
     !!
     !! Interface:
     !!   method1 -> Brief Description of class methods. Very short. Details should
     !!     be contained within comment above procedure definition.
     !!
     type myType
       Type Definition
     end type

#. If type you are writing can be build using a dictionary (usually from user
   input), include additional section with a sample input dictionary. Refer
   to dictionary input syntax :ref:`dictSyntax` ::

      !!
      !! Brief description of the type
      !!
      !! Detailed description of the type
      !!
      !! Public Members:
      !!   A -> Brief description of members. May be missing if there are none
      !!
      !! Private Members:
      !!   B -> Brief Description of members. May be missing if there are none
      !!
      !! Interface:
      !!   method1 -> Brief Description of class methods. Very short. Details should
      !!     be contained within comment above procedure definition.
      !!
      !! Sample Dictionary Input:
      !!   genericName {
      !!     this isMandatoryEntry;
      !!     canAlsoBeNumber 1;
      !!     orList (1 2 3 4 5);
      !!     # Hashes MarkOptionalEntries; #
      !!     # subDict { <name of type that will be build with this subdictionary>} #
      !!   }
      !!
      type myType
        Type Definition
      end type

#. Try to have a look at the code you wrote and just try to make it look pretty.
   Try to go back to your code after a break and try to spot places that seem
   unclear or confusing and improve them.

#. **Comment your Code!** Provide explanations for what given sections are doing,
   explain parts of the algorithms that may be confusing. To see what I mean
   you can have a look at the following code sample. I would argue it is easier
   to understand with comments then without them::

       function sampleLegendre_P1(P1,rand) result(x)
         real(defReal), intent(in)    :: P1
         class(RNG), intent(inout)    :: rand
         real(defReal)                :: x
         real(defReal)                :: P1_loc
         real(defReal)                :: threshold
         integer(shortInt)            :: Low, Top, exec
         integer(shortInt), parameter :: UNIFORM = 1, LIN = 2, DELTA = 3
         character(100), parameter :: Here = 'sampleLegendre_P1 ( legendrePoly_func.f90)'

         ! Make local copy of P1 coeff. Take abs() to simplify code
         ! -ve P1 will be inverted at the end.
         P1_loc = abs(P1)

         ! Depending on whether P1 > 1 determine treshold and associated PDF for the mixing method
         ! For further details refer to Lux and Koblinger APPENDIX 3D
         ! If random number < threshold then Top is used.
         if ( P1_loc < ONE) then
           threshold = P1_loc
           Top = LIN
           Low = UNIFORM

         else if( P1_loc <= 3.0_defReal) then
           threshold = 0.5 * (P1_loc - ONE)
           Top = DELTA
           Low = LIN

         else
           call fatalError(Here,'P1 must have absolute value < 3.0')
           ! Avoid warnings
           threshold = ONE
           Top = 0
           Low = 0

         end if

         ! Use mixing method with the calculated Threshold
         if ( rand % get() < threshold ) then
           exec = Top
         else
           exec = Low
         end if

         ! Sample from UNIFORM ( PDF = 0.5); LIN ( PDF = 0.5 + 0.5 *mu) or DELTA ( PDF = DELTA(mu-1))
         select case(exec)
           case (UNIFORM)
             x = TWO * rand % get() - ONE

           case (LIN)
             ! Need to solve CDF(x) = 0.25 * x^2 + 0.5 * x + 0.25 = (0.5*x+0.5)^2)
             x = TWO * sqrt(rand % get()) - ONE

           case (DELTA)
             x = ONE

           case default
             call fatalError(Here,'This should never happen. WTF?')
             x = ZERO

         end select

         ! Invert result if P1 is -ve
         if ( P1 < ZERO ) x = -x

       end function sampleLegendre_P1
