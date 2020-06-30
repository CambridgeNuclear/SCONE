module surfaceShelf_class

  use numPrecision
  use genericProcedures,   only : fatalError, numToChar
  use dictionary_class,    only : dictionary
  use intMap_class,        only : intMap
  use surface_inter,       only : surface
  use surfaceFactory_func, only : new_surface_ptr

  implicit none
  private

  !!
  !! Small, local container to store polymorphic surfaces in a single array
  !!
  !! Public members:
  !!   name -> Name of the surface
  !!   ptr  -> Pointer to the surface
  !!
  type :: surfaceBox
    character(nameLen)      :: name = ''
    class(surface), pointer :: ptr => null()
  end type surfaceBox

  !!
  !! Storage space for surfaces defined in the geometry
  !!
  !! Sample dictionary input:
  !!   surfaces {
  !!    surf1 { <surface definition> }
  !!    surf2 { <surface definition> }
  !!    surf3 { <surface definition> }
  !!    ...
  !!    }
  !!
  !! Private Members:
  !!   surfaces -> Array to store pointers to polymorphic surfaces
  !!   idMap -> Map between surface ID and corresponding index
  !!
  !! Interface:
  !!   init -> Initialise from a dictionary
  !!   surfPtr -> Return pointer to a surface given its index
  !!   surfIdx -> Return index of a surface given its id
  !!   surfID  -> Return id of a surface given its idx
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Becouse surfaces are stored as pointers, calling `kill` is crutial to prevent
  !!   memory leaks. TODO: Add `final` procedure here ?
  !!
  type, public :: surfaceShelf
    private
    type(surfaceBox), dimension(:), allocatable :: surfaces
    type(intMap)                                :: idMap

  contains
    procedure :: init
    procedure :: surfPtr
    procedure :: surfIdx
    procedure :: surfId
    procedure :: kill
  end type surfaceShelf

contains

  !!
  !! Load surfaces into shelf
  !!
  !! Args:
  !!   dict [in] -> Dictionary with subdictionaries that contain surface definitions
  !!
  !! Errors:
  !!   fatalError if there are clashes in surface ID
  !!
  subroutine init(self, dict)
    class(surfaceShelf), intent(inout)            :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: names
    integer(shortInt)                             :: i, id, idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'init (surfaceShelf_class.f90)'

    ! Get all keys for subdictionaries
    call dict % keys(names, 'dict')

    ! Allocate space
    allocate (self % surfaces(size(names)))

    ! Build surfaces
    do i = 1, size(names)
      self % surfaces(i) % name = names(i)
      self % surfaces(i) % ptr => new_surface_ptr(dict % getDictPtr(names(i)))
      id = self % surfaces(i) % ptr % id()

      ! Add ID to the map detecting any conflicts
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) then
        call fatalError(Here,'Surfaces '//trim(names(i))// ' & '//&
                             trim(self % surfaces(idx) % name)//&
                             ' have the same ID: '//numToChar(id))

      else
        call self % idMap % add(id, i)

      end if
    end do

  end subroutine init

  !!
  !! Return pointer to the surface indicated by index
  !!
  !! Args:
  !!   idx [in] -> Index of the surface
  !!
  !! Result:
  !!   Pointer to a surface under index idx
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a surface (is out-of-bounds)
  !!
  function surfPtr(self, idx) result (ptr)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    class(surface), pointer         :: ptr
    character(100), parameter :: Here = 'surfPtr (surfaceShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % surfaces)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % surfaces)))
    end if

    ! Return pointer
    ptr => self % surfaces(idx) % ptr

  end function surfPtr

  !!
  !! Return IDX of a surface with ID
  !!
  !! Args:
  !!   id [in] -> Id of the surface
  !!
  !! Result:
  !!   Index of a surface with ID
  !!
  !! Error:
  !!   fatalError if there is not surface with ID
  !!
  function surfIdx(self, id) result(idx)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: id
    integer(shortInt)               :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'surfIdx (surfaceShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)

    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'There is no surface with ID: '//numToChar(id))
    end if

  end function surfIdx

  !!
  !! Return ID of the surface with index
  !!
  !! Args:
  !!   idx [in] -> Index of the surface
  !!
  !! Result:
  !!   ID of the surface under index
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a surface (is out-of-bounds)
  !!
  function surfId(self, idx) result(id)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: id
    character(100), parameter :: Here = 'surfID (surfaceShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % surfaces)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % surfaces)))
    end if

    id = self % surfaces(idx) % ptr % id()

  end function surfId

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(surfaceShelf), intent(inout) :: self
    integer(shortInt)                  :: i

    if (allocated(self % surfaces)) then
      do i = 1, size(self % surfaces)
        call self % surfaces(i) % ptr % kill()
      end do

      deallocate(self % surfaces)
    end if
    
    call self % idMap % kill()

  end subroutine kill

end module surfaceShelf_class
