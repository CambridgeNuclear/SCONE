!!
!! Module to build nuclear data classes and store information (names + pointers) about them
!!
module nuclearDataRegistry_mod

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Abstract interfaces
  use nuclearData_inter,              only : nuclearData
  use transportNuclearData_inter,     only : transportNuclearData
  use perNuclideNuclearDataCE_inter,  only : perNuclideNuclearDataCE
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Individual implementations
  use byNucNoMT_class,   only : byNucNoMT
  use byNucMT_class,     only : byNucMT
  use isotropicMG_class, only : isotropicMG
  use transMG_class,     only : transMG
  use P1MG_class,        only : P1MG

  implicit none
!  private

  ! *** ADD NAME OF A NEW NUCLEAR DATA HERE ***!
  ! List that contains all accaptable types of nuclear data
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_nuclearData = [ 'byNucNoMT  ', &
                                                                        'byNucMT    ', &
                                                                        'isotropicMG', &
                                                                        'transMG    ', &
                                                                        'P1MG       ']
  !!
  !! Helper structure to store pointers to multiple nuclear data
  !! objects in an array
  !!
  type, private :: nuclearDataBox
    character(nameLen)          :: name
    class(nuclearData), pointer :: data  => null()
  end type nuclearDataBox


  !public :: new_nuclearData_ptr

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Module Components
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  type(nuclearDataBox),dimension(:),allocatable :: nucData


contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Builds and allocates an instance of nuclear data from dictionary
  !!
  function new_nuclearData_ptr(dict, type, matNames) result(new)
    class(dictionary), intent(in)              :: dict
    character(nameLen),intent(in)              :: type
    character(nameLen),dimension(:),intent(in) :: matNames
    class(nuclearData), pointer                :: new
    character(100),parameter :: Here = 'new_nuclearData_ptr (nuclearDataFactory_func.f90)'

    ! Allocate approperiate subclass of nuclear data
    ! *** ADD CASE STATEMENT FOR A NEW NUCLEAR DATA BELOW ***!
    select case(type)
      case('byNucNoMT')
        ! Allocate and initialise
        allocate( byNucNoMT :: new)

      case('byNucMT')
        ! Allocate and initialise
        allocate( byNucMT :: new)

      case('isotropicMG')
        ! Allocate and initialise
        allocate( isotropicMG :: new)

      case('transMG')
        ! Allocate and initialise
        allocate( transMG :: new)

      case('P1MG')
        ! Allocate and initialise
        allocate( P1MG :: new)

      case default
        print *, AVALIBLE_nuclearData
        call fatalError(Here, 'Unrecognised type of nuclearData: ' // trim(type))

    end select

    ! Initialise an instance of data
    call new % init(dict, matNames)

  end function new_nuclearData_ptr

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Public Interface of the module
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Given a dictionary this subroutine builds all nuclear data representations
  !!
  subroutine nuclearData_buildMaterials(dict)
    class(dictionary), intent(in) :: dict
    class(dictionary), pointer    :: handlesDict   => null()
    class(dictionary), pointer    :: materialsDict => null()

    ! Obtain dictionaries
    handlesDict   => dict % getDictPtr('handles')
    materialsDict => dict % getDictPtr('materials')

    ! Loop over all handles and build data representations and handle name -> handleIdx map

    ! Loop over materials and load them into charMap (be carefull to preserve order).

    ! Verify material order for all data types

  end subroutine nuclearData_buildMaterials



end module nuclearDataRegistry_mod
