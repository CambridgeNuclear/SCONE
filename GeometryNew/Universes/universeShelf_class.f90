module universeShelf_class

  use numPrecision
  use genericProcedures,    only : fatalError, numToChar
  use dictionary_class,     only : dictionary
  use intMap_class,         only : intMap
  use charMap_class,        only : charMap
  use surfaceShelf_class,   only : surfaceShelf
  use cellShelf_class,      only : cellShelf

  ! Universe objects
  use universe_inter,       only : universe
  use universeFactory_func, only : new_universe_ptr
  use uniFills_class,       only : uniFills

  implicit none
  private

  !!
  !! Small, local container to store polymorphic universes in an array
  !!
  !! Public members:
  !!   name -> Name of the universe
  !!   ptr  -> Pointer to the universe
  !!
  type :: uniBox
    character(nameLen)       :: name = ''
    class(universe), pointer :: ptr => null()
  end type uniBox

  !!
  !! Storage space for universes defined in the geometry
  !!
  !! Sample Dictionary Input:
  !!   universes {
  !!     uni1 { <universe definition>}
  !!     uni2 { <universe definition>}
  !!     ...
  !!   }
  !!
  !! Private Members:
  !!   unis -> Array with pointers to different universes
  !!   idMap -> Map between uniId and uniIdx
  !!
  !! Interface:
  !!   init    -> Initialise and build uniFills
  !!   getPtr  -> Get pointer to a universe given by its index
  !!   getIdx  -> Get uniIdx of a universe given by uniId
  !!   getId   -> Get uniId of a universe given by uniIdx
  !!   getSize -> Return the number of universes (max uniIdx)
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Becoue universes are stored as pointers, calling `kill` is crucial
  !!   to prevent memory leaks. TODO: Add `final` procedure here?
  !!
  type, public :: universeShelf
    private
    type(uniBox), dimension(:), allocatable :: unis
    type(intMap)                            :: idMap
  contains
    procedure :: init
    procedure :: getPtr
    procedure :: getIdx
    procedure :: getId
    procedure :: getSize
    procedure :: kill
  end type universeShelf

contains

  !!
  !! Initialise universeShelf
  !!
  !! Builds all definitions and load fill information into uniFills
  !!
  !! Args:
  !!   fills [out]   -> `uniFills` that contains fill information for every universe
  !!   dict [in]     -> Dictionary with universe definitions
  !!   cells [inout] -> `cellShelf` with used defined cells
  !!   surfs [inout] -> `surfaceShelf` with user defined surfaces
  !!   mats [in]     -> Map of material names to matIdx
  !!
  !! Errors:
  !!   fatalError if there are multiple universes with the same id
  !!
  subroutine init(self, fills, dict, cells, surfs, mats)
    class(universeShelf), intent(inout) :: self
    type(uniFills), intent(out)         :: fills
    class(dictionary), intent(in)       :: dict
    type(cellShelf), intent(inout)      :: cells
    type(surfaceShelf), intent(inout)   :: surfs
    type(charMap), intent(in)           :: mats
    character(nameLen), dimension(:), allocatable :: names
    integer(shortInt)                             :: i, id, idx
    integer(shortInt), dimension(:), allocatable  :: fillInfo
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter    :: Here = 'init (universeShelf_class.f90)'

    ! Get all universe names
    call dict % keys(names, 'dict')

    ! Allocate space and uniFills size
    allocate (self % unis(size(names)))
    call fills % init(size(names))

    ! Build universes
    do i = 1, size(names)
      ! Build new universe
      self % unis(i) % name = names(i)

      ! Build new interface
      call new_universe_ptr(self % unis(i) % ptr, &
                            fillInfo, &
                            dict % getDictPtr(names(i)),&
                            cells, surfs, mats)

      ! Add ID to map detecting any conflicts
      id = self % unis(i) % ptr % id()
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) then
        call fatalError(Here,'Universes '//trim(names(i))// ' & '//&
                              trim(self % unis(idx) % name)//&
                             ' have the same ID: '//numToChar(id))
      else
        call self % idMap % add(id, i)

      end if

      ! Store content info in fills
      call fills % addUniverse(i, id, fillInfo)

    end do

    ! Finish build
    call fills % finishBuild(self % idMap)

  end subroutine init

  !!
  !! Return pointer to the universe indicated by idx
  !!
  !! Args:
  !!   idx [in] -> Index of the universe
  !!
  !! Result:
  !!   Pointer to the universe under index idx
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a universe (is out-of-bounds)
  !!
  function getPtr(self, idx) result (ptr)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    class(universe), pointer         :: ptr
    character(100), parameter :: Here = 'getPtr (universeShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % unis)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % unis)))
    end if

    ! Return pointer
    ptr => self % unis(idx) % ptr

  end function getPtr

  !!
  !! Return IDX of a universe with the ID
  !!
  !! Args:
  !!   id [in] -> Id of the universe
  !!
  !! Result:
  !!   Index of the universe with the given Id
  !!
  !! Error:
  !!   fatalError if there is no universe with the given Id
  !!
  function getIdx(self, id) result(idx)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: id
    integer(shortInt)                :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'getIdx (universeShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)

    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'There is no universe with ID: '//numToChar(id))
    end if

  end function getIdx

  !!
  !! Return ID of the universe given by index
  !!
  !! Args:
  !!   idx [in] -> Index of the universe
  !!
  !! Result:
  !!   ID of the universe under the index
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a universe (is out-of-bounds)
  !!
  function getId(self, idx) result(id)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    integer(shortInt)                :: id
    character(100), parameter :: Here = 'getId (universeShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % unis)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % unis)))
    end if

    id = self % unis(idx) % ptr % id()

  end function getId

  !!
  !! Return size of the shelf
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of universes on the shelf
  !!
  elemental function getSize(self) result(N)
    class(universeShelf), intent(in) :: self
    integer(shortInt)                :: N

    N = size(self % unis)

  end function getSize

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(universeShelf), intent(inout) :: self

    if (allocated(self % unis)) deallocate(self % unis)
    call self % idMap % kill()

  end subroutine kill

end module universeShelf_class
