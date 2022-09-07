!!
!! Object Registry that contains all defined nuclearDatabases
!!
!! Stores definitions and memory for all defined databases
!! Allows to obtain a pointer to a database
!! Serves as nuclearDatabase factory
!!
!! Every database can exist in three states:
!!   Unmade -> Name and Definition is loaded. class(nuclearDatabase) is unallocated
!!   Made   -> Class(nuclearDatabase) is allocated and initialised
!!   Active -> Database is Made and associated with one of ND TYPES
!!
!! Available ND TYPES:
!!   CE_NEUTRON
!!   MG_NEUTRON
!!   MG_PHOTON
!!
!! Private members:
!!   databases           -> Array with defined databases (name, definition,
!!   databaseNameMap     -> CharMap that maps database name to index in "databases"
!!   active_ceNeutron    -> Pointer to active CE Neutron database
!!   activeIdx_ceNeutron -> Index of active CE Neutron database in "databases"
!!   active_mgNeutron    -> Pointer to active MG Neutron database
!!   activeIdx_mgNeutron -> Index of active MG Neutron database in "databases"
!!
!! Interface:
!!   init     -> Initialise Nuclear Database Registry
!!   make     -> Use a definition of a ND to allocate it and load the data
!!   clean    -> Purge all the data and deallocate a ND
!!   activate -> Associate a ND with given type of cache and load active materials
!!   display  -> Display activated and defined NDs to console
!!   kill     -> Return to uninitialised state
!!   getNeutronCE -> Get pointer to the active neutron CE Nuclear Database
!!   getNeutronMG -> Get pointer to the active neutron MG Nuclear Database
!!   get          -> Return pointer to active Nuclear Data given particle type or name
!!   getMatNames  -> Returns pointer to charMap of material names and matIdx from materialMenu
!!
!! Note:
!!   To add new nuclearDatabase:
!!     1) Add new "use.." statment
!!     2) Add name of the database to AVAILABLE_NUCLEAR_DATABASES
!!     3) Add new entry in "select case" in "new_nuclearDatabase" procedure
!!     4) Use full name of the database e.g. "aceNeutronDatabase" instead of "aceNeutron"
!!
!!   To add new ND TYPE
!!     1) Create New entries for activeIdx_ and active_
!!     2) In "display" add new entry in "ACTIVE DATABASES" and Add new if(...) cycle to
!!        loop that lits unused databases
!!     3) Add new entry at the end of "kill" subroutine
!!     4) Define new Parameter for the data e.g. "CE_NEUTRON"
!!     5) Add new entry in "select case" in "activate" procedure
!!     6) Add new entry to "Available ND TYPES:" in this doc comment
!!
module nuclearDataReg_mod

  use numPrecision
  use universalVariables,    only : P_NEUTRON_CE, P_NEUTRON_MG, P_PHOTON_MG, P_MATERIAL_MG
  use genericProcedures,     only : fatalError, numToChar, printParticleType
  use charMap_class,         only : charMap
  use dictionary_class,      only : dictionary

  ! Nuclear Data Interfaces & Classes
  use nuclearDatabase_inter,   only : nuclearDatabase
  use ceNeutronDatabase_inter, only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use mgNeutronDatabase_inter, only : mgNeutronDatabase, mgNeutronDatabase_CptrCast
  use mgIMCDatabase_inter,     only : mgIMCDatabase, mgIMCDatabase_CptrCast
  use materialMenu_mod,        only : mm_init => init, mm_kill => kill, mm_nMat => nMat,&
                                      mm_nameMap => nameMap

  ! Implemented Nuclear Databases
  ! Neutron CE
  use aceNeutronDatabase_class,       only : aceNeutronDatabase
  use aceNeutronDatabaseUni_class,    only : aceNeutronDatabaseUni
  use aceNeutronDatabaseUniIdx_class, only : aceNeutronDatabaseUniIdx

  ! Neutron MG
  use baseMgNeutronDatabase_class, only : baseMgNeutronDatabase

  ! Photon MG
  use baseMgIMCDatabase_class,     only : baseMgIMCDatabase

  implicit none
  private

  !!
  !! Local helper container to store polymorphic instances on nuclearDatabases
  !!
  !! Public Members:
  !!   nd -> Polymorphic Nuclear Database
  !!
  type, private :: ndBox
    character(nameLen)                  :: name
    type(dictionary)                    :: def
    class(nuclearDatabase), allocatable :: nd
  end type


  !! Public Interface
  public :: init
  public :: make
  public :: clean
  public :: activate
  public :: display
  public :: kill
  public :: getNeutronCE
  public :: getNeutronMG
  public :: getIMCMG
  public :: get
  public :: getMatNames

  !! Procedures grouped under name "get"
  interface get
    module procedure :: get_byType
    module procedure :: get_byName
  end interface


  !! Parameters
  character(nameLen), dimension(*), parameter :: AVAILABLE_NUCLEAR_DATABASES = &
                                                ['aceNeutronDatabase      ', &
                                                 'baseMgNeutronDatabase   ', &
                                                 'baseMgIMCDatabase       ', &
                                                 'aceNeutronDatabaseUni   ', &
                                                 'aceNeutronDatabaseUniIdx']

  !! Members
  type(ndBox),dimension(:),allocatable,target :: databases
  type(charMap)                               :: databaseNameMap

  class(ceNeutronDatabase), pointer :: active_ceNeutron => null()
  integer(shortInt)                 :: activeIdx_ceNeutron = 0

  class(mgNeutronDatabase), pointer :: active_mgNeutron => null()
  integer(shortInt)                 :: activeIdx_mgNeutron = 0

  class(mgIMCDatabase),     pointer :: active_mgIMC     => null()
  integer(shortInt)                 :: activeIdx_mgIMC     = 0

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
  !!     // Keyword is ND name, Contents of nested dictionary its settings
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
    class(dictionary), intent(in)                 :: dict
    class(dictionary), pointer                    :: handles
    character(nameLen), dimension(:), allocatable :: dataNames
    integer(shortInt)                             :: i

    ! Get pointer to database handles
    handles => dict % getDictPtr('handles')

    ! Get names of databases
    call handles % keys(dataNames, 'dict')

    ! Allocate space
    allocate(databases(size(dataNames)))

    ! Load definitions
    ! Associate names with idx's in Map
    do i=1,size(databases)
      databases(i) % name = dataNames(i)
      databases(i) % def  = handles % getDictPtr(dataNames(i)) ! Note deep copy
      call databaseNameMap % add(dataNames(i), i)
    end do

    ! Load Materials
    call mm_init(dict % getDictPtr('materials'))

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
  !!   If database if already made it has NO EFFECT
  !!
  subroutine make(name, silent)
    character(nameLen), intent(in)         :: name
    logical(defBool), optional, intent(in) :: silent
    logical(defBool)                       :: silent_loc
    integer(shortInt)                      :: idx
    character(nameLen)                     :: type
    class(nuclearDatabase),pointer         :: ptr
    character(100),parameter :: Here = 'make (nuclearDataReg_mod.f90)'

    ! Process optional arguments
    if(present(silent)) then
      silent_loc = silent
    else
      silent_loc = .false.
    end if

    ! Get index
    idx = databaseNameMap % getOrDefault(name, 0)
    if(idx == 0 ) then
      call fatalError(Here, trim(name)//' is was not defined. Cannot make it!')
    else if(idx < 0) then
      call fatalError(Here, '-ve idx from databaseNameMap. Quite immpossible. WTF?')
    end if

    ! Quit if already has been allocated
    if(allocated(databases(idx) % nd)) return

    ! Build Nuclear Database
    call databases(idx) % def % get(type, 'type')
    call new_nuclearDatabase(databases(idx) % nd, type)

    ! Initialise
    ptr => databases(idx) % nd
    call databases(idx) % nd % init( databases(idx) % def, ptr, silent = silent_loc)

  end subroutine make

  !!
  !! Deallocate the selected database
  !!
  !! Args:
  !!   name [in] -> Name of a Nuclear Database to clean
  !!
  !! Errors:
  !!   NO EFFECT if Database under name is not present
  !!   NO EFFECT if Database is not made
  !!
  subroutine clean(name)
    character(nameLen), intent(in) :: name
    integer(shortInt)              :: idx

    idx = databaseNameMap % getOrDefault(name, 0)
    if (idx < 1) return

    if(allocated(databases(idx) % nd)) then
      call databases(idx) % nd % kill()
      deallocate(databases(idx) % nd)
    end if

  end subroutine clean

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
  !!   fatalError if name is not present
  !!
  !!
  subroutine activate(type, name, activeMat, silent)
    integer(shortInt), intent(in)               :: type
    character(nameLen), intent(in)              :: name
    integer(shortInt), dimension(:) ,intent(in) :: activeMat
    logical(defBool), optional, intent(in)      :: silent
    logical(defBool)                            :: silent_loc
    integer(shortInt)                           :: idx
    class(nuclearDatabase), pointer             :: ptr
    character(100), parameter :: Here = 'activate (nuclearDataReg_mod.f90)'

    ! Process Optional Arguments
    silent_loc = .false.
    if (present(silent)) silent_loc = silent

    ! Get index
    idx = databaseNameMap % getOrDefault(name, 0)
    if(idx == 0 ) then
      call fatalError(Here, trim(name)//' is was not defined. Cannot activate it!')
    else if(idx < 0) then
      call fatalError(Here, '-ve idx from databaseNameMap. Quite immpossible. WTF?')
    end if

    ! Make if it is not already made
    if(.not.allocated(databases(idx) % nd)) call make(name, silent = silent_loc)

    ! Activate
    call databases(idx) % nd % activate(activeMat)
    ptr => databases(idx) % nd

    ! Register as active
    ! This is a bit of a messy code. Could be better. Blame me. - MAK
    select case(type)
      case(P_NEUTRON_CE)
        activeIdx_ceNeutron = idx
        active_ceNeutron => ceNeutronDatabase_CptrCast(ptr)
        if(.not.associated(active_ceNeutron)) then
          call fatalError(Here,trim(name)//' is not database for CE neutrons')
        end if

      case(P_NEUTRON_MG)
        activeIdx_mgNeutron = idx
        active_mgNeutron => mgNeutronDatabase_CptrCast(ptr)
        if(.not.associated(active_mgNeutron)) then
          call fatalError(Here,trim(name)//' is not database for MG neutrons')
        end if

      case(P_PHOTON_MG)
        activeIdx_mgIMC = idx
        active_mgIMC => mgIMCDatabase_CptrCast(ptr)
        if(.not.associated(active_mgIMC)) then
          call fatalError(Here,trim(name)//' is not database for MG IMC')
        end if

      case default
        call fatalError(Here,'Unrecognised type of data to activate. Check parameters. Got: '//&
                              numToChar(type))
    end select

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
    integer(shortInt)  :: idx
    character(nameLen) :: activeName

    print '(A)',repeat('/\',30)
    print '(A)', "ACTIVE DATABASES:"

    ! CE NEUTRON
    activeName = 'NONE'
    idx = activeIdx_ceNeutron
    if(idx /= 0) activeName = databases(idx) % name
    print '(A)', "  CE NEUTRON DATA: " // trim(activeName)

    ! MG NEUTRON
    activeName = 'NONE'
    idx = activeIdx_mgNeutron
    if(idx /= 0) activeName = databases(idx) % name
    print '(A)', "  MG NEUTRON DATA: " // trim(activeName)

    ! MG IMC
    activename = 'NONE'
    idx = activeIdx_mgIMC
    if(idx /= 0) activeName = databases(idx) % name
    print '(A)', "  MG IMC DATA: "     // trim(activeName)

    ! INACTIVE DATABASES
    print '(A)', "INACTIVE DATABASES:"
    do idx=1,size(databases)
      if(idx == activeIdx_mgNeutron) cycle
      if(idx == activeIdx_ceNeutron) cycle
      if(idx == activeIdx_mgIMC)     cycle

    end do
    print '(A)',repeat('\/',30)
  end subroutine display

  !!
  !! Return to uninitialised state
  !!
  subroutine kill()
    integer(shortInt) :: it
    !! Clean all databases
    it = databaseNameMap % begin()
    do while (it /= databaseNameMap % end())
      call clean(databaseNameMap % atKey(it))
      it = databaseNameMap % next(it)

    end do

    !! Take care of databases array
    if(allocated(databases)) then
      do it =1,size(databases)
        call databases(it) % def % kill()
      end do
      deallocate(databases)
    end if

    !! Return pointers to active databases to initial state
    ! CE NEUTRON
    activeIdx_ceNeutron = 0
    active_ceNeutron => null()

    ! MG NEUTRON
    activeIdx_mgNeutron = 0
    active_mgNeutron => null()

    ! MG IMC
    activeIdx_mgIMC     = 0
    active_mgIMC     => null()

  end subroutine kill

  !!
  !! Return pointer to an active Neutron CE Database
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   ceNeutronDatabase class pointer
  !!
  !! Errors:
  !!   If there is no active database returns NULL ptr
  !!
  function getNeutronCE() result(ptr)
    class(ceNeutronDatabase), pointer :: ptr

    ptr => active_ceNeutron

  end function getNeutronCE

  !!
  !! Return pointer to an active Neutron MG Database
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   mgNeutronDatabase class pointer
  !!
  !! Errors:
  !!   If there is no active database returns NULL ptr
  !!
  function getNeutronMG() result(ptr)
    class(mgNeutronDatabase), pointer :: ptr

    ptr => active_mgNeutron

  end function getNeutronMG

  !!
  !! Return pointer to an active IMC MG Database
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   mgIMCDatabase class pointer
  !!
  !! Errors:
  !!   If there is no active database returns NULL ptr
  !!
  function getIMCMG() result(ptr)
    class(mgIMCDatabase), pointer :: ptr

    ptr => active_mgIMC

  end function getIMCMG

  !!
  !! Return pointer to an active Nuclear Database given particle type
  !!
  !! Args:
  !!   type [in]  -> Particle type
  !!   where [in] -> Optional, Location of error message
  !!
  !! Result:
  !!   nuclearDatabaseclass pointer
  !!
  !! Errors:
  !!   fatalError if there no activa database or type is invalid
  !!
  function get_byType(type, where) result(ptr)
    integer(shortInt), intent(in)     :: type
    class(nuclearDatabase), pointer   :: ptr
    character(*),optional, intent(in) :: where
    character(100), parameter         :: Here = 'get_byType (nuclearDataReg_mod.f90)'

    select case(type)
      case(P_NEUTRON_CE)
        ptr => getNeutronCE()

      case(P_NEUTRON_MG)
        ptr => getNeutronMG()

      case(P_PHOTON_MG)
        ptr => getIMCMG()

      case(P_MATERIAL_MG)
        ! Currently only used for ISMC so point to same database
        ptr => getIMCMG()

      case default
        ptr => null()
    end select

    ! Throw error if somthing went wrong
    if(.not.associated(ptr) .and. present(where)) then
      call fatalError(Where, "There is no data for particle: "//printParticleType(type))

    else if(.not.associated(ptr)) then
      call fatalError(Here, "There is no data for particle: "//printParticleType(type))

    end if

  end function get_byType

  !!
  !! Return pointer to a Nuclear Databse given its name
  !!
  !! Args:
  !!   name [in]  -> Name of the database
  !!   where [in] -> Optional, Location of error message
  !!
  !! Result:
  !!   nuclearDatabaseclass pointer
  !!
  !! Errors:
  !!   fatalError if name is invalid
  !!
  function get_byName(name, where) result(ptr)
    character(*), intent(in)         :: name
    class(nuclearDatabase), pointer  :: ptr
    character(*),optional,intent(in) :: where
    character(nameLen)               :: name_loc
    integer(shortInt)                :: idx
    character(100), parameter        :: Here = 'get_byType (nuclearDataReg_mod.f90)'

    name_loc = name
    idx = databaseNameMap % getOrDefault(name_loc, -1)

    if(idx == -1 .and. present(where)) then
      call fatalError(where, name // " was not found among databases")
      ptr => null() ! Avoid warning

    else if (idx == -1) then
      call fatalError(Here, name // " was not found among databases")
      ptr => null() ! Avoid warning

    else
      ptr => databases(idx) % nd
    end if

  end function get_byName



  !!
  !! Return pointer to charMap of materialNames to matIdx from MaterialMenu
  !!
  !! It existis to hide the existance of materialMenu outside of nucleardata
  !! and limit the interface to fewer modules
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Pointer to charMap with materialNames and matIdx
  !!
  !! Errors:
  !!   fatalError if Material Menu was not yet initialised
  !!
  function getMatNames() result(ptr)
    type(charMap), pointer :: ptr
    character(100),parameter :: Here = 'getMatNames (nuclearDataReg_mod.f90)'

    if (mm_nMat() == 0) call fatalError(Here, "Material Definitions are empty. Has Nuclear Data been initialised?")
    ptr => mm_nameMap

  end function getMatNames


  !!
  !! Allocates the database to a specified type
  !!
  !! This is Factory procedure for nuclearDatabases
  !!
  !! Args:
  !!   database [inout] -> allocatable class(nuclearDatabase) to be allocated
  !!   type [in]        -> character that specifies type to be allocated
  !!
  !! Errors:
  !!   fatalError if type is not recognised. AVALIABLE_NUCLEAR_DATABASES will also be printed
  !!   If database is allocated on entry it will be killed and reallocated
  !!
  subroutine new_nuclearDatabase(database, type)
    class(nuclearDatabase), allocatable , intent(inout) :: database
    character(nameLen), intent(in)                      :: type
    integer(shortInt)                                   :: i
    character(100), parameter :: Here = 'new_nuclearDatabase (nuclearDataReg_mod.f90)'

    ! Kill if needed
    if(allocated(database)) then
      call database % kill()
      deallocate(database)
    end if

    ! Allocate to required type
    select case(type)
      case('aceNeutronDatabase')
        allocate(aceNeutronDatabase :: database)

      case('baseMgNeutronDatabase')
        allocate(baseMgNeutronDatabase :: database)

      case('baseMgIMCDatabase')
        allocate(baseMgIMCDatabase :: database)

      case('aceNeutronDatabaseUni')
        allocate(aceNeutronDatabaseUni :: database)

      case('aceNeutronDatabaseUniIdx')
        allocate(aceNeutronDatabaseUniIdx :: database)

      case default
        ! Print available nuclear database types
        print '(A)', "<><><><><><><><><><><><><><><><><><><><>"
        print '(A)', "Available Nuclear Databases:"
        do i=1,size(AVAILABLE_NUCLEAR_DATABASES)
          print '(A)', AVAILABLE_NUCLEAR_DATABASES(i)
        end do

        ! Throw FatalError
        call fatalError(Here,trim(type) //' is not valid nuclearDatabase. See list Above')
    end select

  end subroutine new_nuclearDatabase


end module nuclearDataReg_mod
