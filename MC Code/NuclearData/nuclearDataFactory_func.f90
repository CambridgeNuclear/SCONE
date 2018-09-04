!!
!! Module to build nuclear data classes and store information (names + pointers) about them
!!
module nuclearDataFactory_func

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
  use transMG_class,      only : transMG

  implicit none
  private

  ! *** ADD NAME OF A NEW NUCLEAR DATA HERE ***!
  ! List that contains all accaptable types of nuclear data
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_nuclearData = [ 'byNucNoMT  ', &
                                                                        'byNucMT    ', &
                                                                        'isotropicMG', &
                                                                        'transMG    ']

  public :: new_nuclearData_ptr

contains

  !!
  !! Builds and allocates an instance of nuclear data from dictionary
  !!
  function new_nuclearData_ptr(dict) result(new)
    class(dictionary), intent(in)         :: dict
    class(nuclearData), pointer           :: new
    character(nameLen)                    :: type
    character(100),parameter :: Here = 'new_nuclearData_ptr (nuclearDataFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of nuclear data
    ! *** ADD CASE STATEMENT FOR A NEW NUCLEAR DATA BELOW ***!
    ! **** AT THE MOMENT ALLOCATE + SELECT TYPE + INIT is very unelegant implementation
    ! **** Will have to be improved
    select case(type)
      case('byNucNoMT')
        ! Allocate and initialise
        allocate( byNucNoMT :: new)
        select type(new)
          type is (byNucNoMT)
            call new % init(dict)
        end select

      case('byNucMT')
        ! Allocate and initialise
        allocate( byNucMT :: new)
        select type(new)
          type is (byNucMT)
            call new % init(dict)
        end select

      case('isotropicMG')
        ! Allocate and initialise
        allocate( isotropicMG :: new)
        select type(new)
          type is (isotropicMG)
            call new % init(dict)
        end select

      case('transMG')
        ! Allocate and initialise
        allocate( transMG :: new)
        select type(new)
          type is (transMG)
            call new % init(dict)
        end select

      case default
        print *, AVALIBLE_nuclearData
        call fatalError(Here, 'Unrecognised type of nuclearData: ' // trim(type))

    end select

  end function new_nuclearData_ptr

end module nuclearDataFactory_func
