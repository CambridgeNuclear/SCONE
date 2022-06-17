module universeFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf

  ! Universe interface
  use universe_inter, only : universe

  ! Universes
  use rootUniverse_class,     only : rootUniverse
  use cellUniverse_class,     only : cellUniverse
  use pinUniverse_class,      only : pinUniverse
  use azimPinUniverse_class,  only : azimPinUniverse
  use latUniverse_class,      only : latUniverse
  implicit none
  private

  ! ** ADD NAME OF NEW UNIVERSE TO THE LIST
  ! List contains acceptable types of universe
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_UNI = ['rootUniverse    ',&
                                                                  'cellUniverse    ',&
                                                                  'pinUniverse     ',&
                                                                  'azimPinUniverse ',&
                                                                  'latUniverse     ']

  ! Public Interface
  public :: new_universe_ptr
  public :: new_universe

contains

  !!
  !! Point a pointer to a new instance of an allocated universe
  !!
  !! Args:
  !!  ptr [out]     -> Pointer to the new universe
  !!  fill [out]    -> Allocatable integer array with filling of diffrent uniqueID
  !!  dict [in]     -> Dictionary with universe definition
  !!  cells [inout] -> Shelf with defined cells
  !!  surfs [inout] -> Shelf  with defined surfaces
  !!  mats [in]     -> Map of material names to matIdx
  !!
  !! Errors:
  !!   fatalError if type of universe is unknown
  !!
  subroutine new_universe_ptr(ptr, fill, dict, cells, surfs, mats)
    class(universe), pointer, intent(out) :: ptr
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    character(nameLen)                                        :: type
    character(100), parameter :: Here = 'new_universe_ptr (universeFactory_func.f90)'

    ! Obtain type of the universe
    call dict % get(type, 'type')

    ! Allocate approperiate universe
    ! ** FOR NEW UNIVERSE ADD CASE STATEMENT HERE ** !
    select case (type)
      case ('rootUniverse')
        allocate(rootUniverse :: ptr)

      case ('cellUniverse')
        allocate(cellUniverse :: ptr)

      case ('pinUniverse')
        allocate(pinUniverse :: ptr)

      case ('azimPinUniverse')
        allocate(azimPinUniverse :: ptr)

      case ('latUniverse')
        allocate(latUniverse :: ptr)

      case default
        print '(A)', 'AVAILABLE UNIVERSES: '
        print '(A)', AVAILABLE_UNI
        call fatalError(Here, 'Unrecognised type of an universe: '//trim(type))

    end select

    ! Initialise universe
    call ptr % init(fill, dict, cells, surfs, mats )

  end subroutine new_universe_ptr

  !!
  !! Allocte an allocatable universe
  !!
  !! Args:
  !!  new [out]     -> Universe to be allocated
  !!  fill [out]    -> Allocatable integer array with filling of diffrent uniqueID
  !!  dict [in]     -> Dictionary with universe definition
  !!  cells [inout] -> Shelf with defined cells
  !!  surfs [inout] -> Shelf  with defined surfaces
  !!  mats [in]     -> Map of material names to matIdx
  !!
  !! Errors:
  !!   fatalError if type of universe is unknown
  !!
  subroutine new_universe(new, fill, dict, cells, surfs, mats)
    class(universe), allocatable, intent(out)                 :: new
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    class(universe), pointer                                  :: temp

    ! Allocate temporary
    call new_universe_ptr(temp, fill, dict, cells, surfs, mats)
    allocate(new, source = temp)

    ! Clean up
    call temp % kill()
    deallocate(temp)

  end subroutine new_universe

end module universeFactory_func
