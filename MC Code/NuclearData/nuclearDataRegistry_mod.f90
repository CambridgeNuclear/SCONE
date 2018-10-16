!!
!! Module to build nuclear data classes and store information (names + pointers) about them
!!
module nuclearDataRegistry_mod

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use charMap_class,     only : charMap

  ! Abstract interfaces
  use nuclearData_inter,              only : nuclearData
  use transportNuclearData_inter,     only : transportNuclearData
  use perNuclideNuclearDataCE_inter,  only : perNuclideNuclearDataCE
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Individual implementations
  use datalessMaterials_class, only : datalessMaterials
  use byNucNoMT_class,         only : byNucNoMT
  use byNucMT_class,           only : byNucMT
  use isotropicMG_class,       only : isotropicMG
  use transMG_class,           only : transMG
  use P1MG_class,              only : P1MG

  implicit none
  private

  ! Module public interface
  public :: build_nuclearData
  public :: kill_nuclearData
  public :: getMatIdx
  public :: getHandlePtr

  ! *** ADD NAME OF A NEW NUCLEAR DATA HERE ***!
  ! List that contains all accaptable types of nuclear data
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_nuclearData = [ 'byNucNoMT        ', &
                                                                        'byNucMT          ', &
                                                                        'isotropicMG      ', &
                                                                        'transMG          ', &
                                                                        'P1MG             ', &
                                                                        'datalessMaterials']
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
  type(charMap)                                 :: matNames
  type(charMap)                                 :: handleNames

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

      case('datalessMaterials')
        ! Allocate
        allocate( datalessMaterials :: new)

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
  subroutine build_NuclearData(dict)
    class(dictionary), intent(in)               :: dict
    class(dictionary), pointer                  :: handlesDict   => null()
    class(dictionary), pointer                  :: materialsDict => null()
    character(nameLen),dimension(:),allocatable :: names
    character(nameLen)                          :: type
    integer(shortInt)                           :: i,j
    character(100),parameter :: Here ='nuclearData_buildMaterials (nuclearDataRegisty_mod.f90)'

    ! Obtain dictionaries
    handlesDict   => dict % getDictPtr('handles')
    materialsDict => dict % getDictPtr('materials')

    ! Retrieve handles names and allocate storage for pointers
    call handlesDict % keys(names)
    allocate(nucData( size(names) ))
    nucData % name = names

    ! Load handles names into charMap
    call handleNames % init(size(names))
    do i=1,size(names)
      call handleNames % add(names(i), i)

    end do

    ! Loop over materials and load them into charMap (be carefull to preserve order).
    call materialsDict % keysDict(names)
    call matNames % init(size(names))
    do i=1,size(names)
      call matNames % add(names(i), i)

    end do

    ! Build nuclear data representations
    do i=1,size(nucData)
      call handlesDict % get(type, nucData(i) % name)
      nucData(i) % data  => new_nuclearData_ptr(materialsDict, type, names)

    end do

    ! Verify material order for all data types
    do i=1,size(nucData)
      do j=1,size(names)
        if( nucData(i) % data % getIdx(names(j)) /= j) then
          call fatalError(Here,'Inconsistent material names for representation: '// &
                                nucData(i) % name)
        end if
      end do
    end do

  end subroutine build_NuclearData

  !!
  !! Given name of the material return its index
  !!
  function getMatIdx(matName) result(idx)
    character(nameLen), intent(in) :: matName
    integer(shortInt)              :: idx
    integer(shortInt), parameter   :: NOT_FOUND = -huge(idx)
    character(100),parameter :: Here ='getMatIdx (nuclearDataRegistry_mod.f90)'

    idx = matNames % getOrDefault(matName, NOT_FOUND)
    if (idx == NOT_FOUND) then
      call fatalError(Here,'Material '// trim(matName) // ' is not defined')

    else if (idx <= 0) then
      call fatalError(Here,'material index is 0 or -ve. WTF?')

    end if

  end function getMatIdx

  !!
  !! Return pointer to nuclear data handle with the given name
  !!
  function getHandlePtr(handleName) result(ptr)
    character(nameLen), intent(in) :: handleName
    class(nuclearData),pointer     :: ptr
    integer(shortInt)              :: idx
    integer(shortInt), parameter   :: NOT_FOUND = -huge(idx)
    character(100),parameter :: Here ='getHandlePtr (nuclearDataRegistry_mod.f90)'

    idx = handleNames % getOrDefault(handleName, NOT_FOUND)
    if (idx == NOT_FOUND) then
      call fatalError(Here,'Handle '// trim(handleName) // ' is not defined')

    else if (idx <= 0) then
      call fatalError(Here,'handle index is 0 or -ve. WTF?')

    end if

    ptr => nucData(idx) % data

  end function getHandlePtr

  !!
  !! Deallocate memory taken by nuclear data
  !!
  subroutine kill_nuclearData()
    integer(shortInt) :: i

    if(allocated(nucData)) then
      ! Deallocate data
      do i=1,size(nucData)
        call nucData(i) % data % kill()
      end do
      ! Deallocate containers
      deallocate(nucData)

    end if

    call matNames % kill()
    call handleNames % kill()

  end subroutine kill_nuclearData

end module nuclearDataRegistry_mod
