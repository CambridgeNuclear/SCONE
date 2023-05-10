!!
!! Module to generate physicsPackages
!!
module physicsPackageFactory_func

  use numPrecision
  use genericProcedures,               only : fatalError
  use dictionary_class,                only : dictionary

  ! Physics Package interface
  use physicsPackage_inter,            only : physicsPackage

  ! Implementations
  use eigenPhysicsPackage_class,       only : eigenPhysicsPackage
  use fixedSourcePhysicsPackage_class, only : fixedSourcePhysicsPackage
  use vizPhysicsPackage_class,         only : vizPhysicsPackage
  use rayVolPhysicsPackage_class,      only : rayVolPhysicsPackage
  use IMCPhysicsPackage_class,         only : IMCPhysicsPackage

  implicit none
  private

  ! *** ADD NAME OF A NEW PHYSICS PACKAGE HERE ***!
  ! List that contains all accaptable types of Physics Packages
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVAILABLE_physicsPackages = [ 'eigenPhysicsPackage      ',&
                                                                             'fixedSourcePhysicsPackage',&
                                                                             'IMCPhysicsPackage        ',&
                                                                             'vizPhysicsPackage        ',&
                                                                             'rayVolPhysicsPackage     ']

  !!
  !! Public interface
  !!
  public :: new_physicsPackage
  public :: new_physicsPackage_ptr

contains

  !!
  !! Return a new allocatable instance of a physicsPackage
  !!
  function new_physicsPackage(dict) result(new)
    class(dictionary), intent(inout)  :: dict
    class(physicsPackage),allocatable :: new
    character(nameLen)                    :: type
    character(100),parameter :: Here = 'new_physicsPackage (physicsPackageFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of physicsPackage
    ! *** ADD CASE STATEMENT FOR A PHYSICS PACKAGE BELOW ***!
    ! **** AT THE MOMENT ALLOCATE + SELECT TYPE + INIT is very unelegant implementation
    ! **** Will have to be improved
    select case(type)
      case('eigenPhysicsPackage')
        ! Allocate and initialise
        allocate( eigenPhysicsPackage :: new)
        select type(new)
          type is (eigenPhysicsPackage)
            call new % init(dict)
        end select

      case('fixedSourcePhysicsPackage')
        ! Allocate and initialise
        allocate( fixedSourcePhysicsPackage :: new)
        select type(new)
          type is (fixedSourcePhysicsPackage)
            call new % init(dict)
        end select

      case('IMCPhysicsPackage')
        ! Allocate and initialise
        allocate( IMCPhysicsPackage :: new)
        select type(new)
          type is (IMCPhysicsPackage)
            call new % init(dict)
        end select


      case('vizPhysicsPackage')
        ! Allocate and initialise
        allocate( vizPhysicsPackage :: new)
        select type(new)
          type is (vizPhysicsPackage)
            call new % init(dict)
        end select

      case('rayVolPhysicsPackage')
        ! Allocate and initialise
        allocate( rayVolPhysicsPackage :: new)
        call new % init(dict)

      case default
        print *, AVAILABLE_physicsPackages
        call fatalError(Here, 'Unrecognised type of Physics Package : ' // trim(type))

    end select

  end function new_physicsPackage

  !!
  !! Returns allocated pointer to Physics Package build with the dictionary dict
  !!
  function new_physicsPackage_ptr(dict) result(new)
    class(dictionary), intent(inout) :: dict
    class(physicsPackage),pointer    :: new
    character(100),parameter  :: Here = 'new_physicsPackage_ptr (physicsPackageFactory_func.f90)'

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_physicsPackage(dict) )

  end function new_physicsPackage_ptr

end module physicsPackageFactory_func
