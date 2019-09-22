!!
!! Object Registry that contains all defined nuclearDatabases
!!
!! Stores definitions and memory for all defined databases
!! Allows to obtain a pointer to a database
!! Swerves as nuclearDatabase factory
!!
!!
!! Private members:
!!
!! Interface:
!!   init     -> Initialise Nuclear Database Registry
!!   make     -> Use a definition of a ND to allocate it and load the data
!!   clean    -> Purge all the data and deallocate a ND
!!   remake   -> Shorthand for clean + make
!!   activate -> Associate a ND with given type of cache and load active materials
!!   display  -> Display activated and defined NDs to console
!!   kill     -> Return to uninitialised state
!!
module nuclearDataReg_mod

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use dictionary_class,      only : dictionary

  ! Nuclear Data Interfaces & Classes
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Local helper container to store polymorphic instances on nuclearDatabases
  !!
  !! Public Members:
  !!   nd -> Polymorphic Nuclear Database
  !!
  type, private :: ndBox
    class(nuclearDatabase), allocatable :: nd
  end type


  !! Public Interface
  public :: init
  public :: make
  public :: clean
  public :: remake
  public :: activate
  public :: display
  public :: kill

  !! Members


contains

  !!
  !! Initialise Nuclear Data Registry
  !!
  !! Copy dictionaries with settings
  !! Initialise materialMenu
  !!
  !! Args:
  !!   dict [in] -> SCONE dictionary with the data
  !!
  !! Errors:
  !!   fatalError if required data is missing in the dictionary
  !!
  !! Sample Dictionary Input:
  !!
  !! nuclearData {
  !!   databases {
  !!     // Kayword is ND name, Contents of nested dictionary its settings
  !!     ce { type aceNeutronDatabase; aceLib /home/SkekSil/MHmmmmmm/data; }
  !!     mg { type basicMgNeutronDatabase; PN P0; }
  !!   }
  !!   materials {
  !!     mat 1 {             // Refer to materialMenu for detail on materialDefinitions
  !!       temp 273;
  !!       xsFile ./Data/xs1;
  !!       composition {
  !!         1001.03 10.46;
  !!         8016.03  5.23;
  !!       }
  !!     }
  !!   }
  !! }
  !!
  !!
  subroutine init(dict)
    class(dictionary), intent(in) :: dict

  end subroutine init

  !!
  !! Allocate and load data into Database indicated by name
  !!
  !! Args:
  !!   name   [in] -> Name of a Nuclear Database to make
  !!   silent [in] -> Optional. Set to .false. to disable console output
  !!
  !! Errors:
  !!   fatalError if Database under name is not present
  !!
  subroutine make(name, silent)
    character(nameLen), intent(in)         :: name
    logical(defBool), optional, intent(in) :: silent
  end subroutine make

  !!
  !! Deallocate the selected database
  !!
  !! Args:
  !!   name [in] -> Name of a Nuclear Database to clean
  !!
  !! Errors:
  !!   NO EFFECT if Database under name is not present
  !!
  subroutine clean(name)
    character(nameLen), intent(in) :: name
  end subroutine clean

  !!
  !! Remakes the selected database
  !!
  !! Shorthand: Calls clean and make consecutively
  !!
  !! See make and clean docs for extra detail
  !!
  subroutine remake(name, silent)
    character(nameLen), intent(in)         :: name
    logical(defBool), optional, intent(in) :: silent
  end subroutine remake

  !!
  !! Associate a database with one of the caches and set active materials
  !!
  !! If database under name in uninitialised it will be made.
  !!
  !! Args:
  !!   type [in]      -> Integer selector of type of cache
  !!   name [in]      -> Name of a database to activate
  !!   activeMat [in] -> Array of active material matIdx
  !!   silent [in]    -> Optional. Set to .false. to disable console output
  !!
  !! Errors:
  !!
  !!
  subroutine activate(type, name, activeMat, silent)
    integer(shortInt), intent(in)               :: type
    character(nameLen), intent(in)              :: name
    integer(shortInt), dimension(:) ,intent(in) :: activeMat
    logical(defBool), optional, intent(in)      :: silent

  end subroutine activate

  !!
  !! Print active and defined databases to console
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  subroutine display()

  end subroutine display

  !!
  !! Return to uninitialised state
  !!
  subroutine kill()

  end subroutine kill
end module nuclearDataReg_mod
