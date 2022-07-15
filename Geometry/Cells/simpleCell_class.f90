module simpleCell_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, hasDuplicates, numToChar
  use dictionary_class,   only : dictionary
  use surfaceShelf_class, only : surfaceShelf
  use surface_inter,      only : surface
  use cell_inter,         only : cell, kill_super => kill

  implicit none
  private

  !!
  !! Local type to hold pointer to a surface together with a halfspace information
  !!
  !! Public Members:
  !!   surfIdx -> Index of the surface (+ve if +ve halspace is used to define the
  !!     cell, -ve otherwise)
  !!   ptr -> Pointer to a surface
  !!
  type :: surfInfo
    integer(shortInt)       :: surfIdx = 0
    class(surface), pointer :: ptr => null()
  end type


  !!
  !! CSG cell defined by intersection of halfspaces
  !!
  !! The simplest CSG cell, that is defined only by the intersection of
  !! surface halfspaces. It contains no union or complement operations
  !!
  !! Sample input dictionary:
  !!   cell {
  !!     type simpleCell;
  !!     id 17;
  !!     surfaces ( -3 4 8); // Each integer is a valid surface ID
  !!     <content info as defined in cellShelf_class >
  !!   }
  !!
  !! Private Members:
  !!   surfaces -> Array that contains, halfspace info, surfIdx and a pointer to the surface.
  !!     Negative halfspace is indicated by surfIdx < 0. Positive by surfIdx > 0
  !!
  !! Interface:
  !!   cell interface
  !!
  type, public, extends(cell) :: simpleCell
    private
    type(surfInfo), dimension(:), allocatable :: surfaces
  contains
    procedure :: init
    procedure :: inside
    procedure :: distance
    procedure :: kill
  end type simpleCell

contains

  !!
  !! Initialise cell
  !!
  !! See cell_inter for details
  !!
  subroutine init(self, dict, surfs)
    class(simpleCell), intent(inout)  :: self
    class(dictionary), intent(in)     :: dict
    type(surfaceShelf), intent(inout) :: surfs
    integer(shortInt), dimension(:), allocatable :: surfIDs
    integer(shortInt)                            :: surfIdx, id, i
    character(100), parameter :: Here = 'init (simpleCell_class.f90)'

    ! Get surface IDs from dictionary
    call dict % get(surfIDs, 'surfaces')

    ! Allocate space
    allocate(self % surfaces(size(surfIDs)))

    ! Load indexes & Pointers
    do i = 1, size(self % surfaces)
      surfIdx = surfs % getIdx( abs(surfIDs(i)))
      self % surfaces(i) % ptr => surfs % getPtr(surfIdx)
      self % surfaces(i) % surfIdx = sign(surfIdx, surfIDs(i))
    end do

    ! Get cell ID
    call dict % get(id, 'id')
    call self % setId(id)

    ! Check surfaces for duplicates
    if (hasDuplicates(abs(self % surfaces % surfIdx))) then
      call fatalError(Here, 'There are repeated surfaces in definition of cell: '//numToChar(id))
    end if

  end subroutine init

  !!
  !! Return TRUE if position is inside the cell
  !!
  !! See cell_inter for details
  !!
  pure function inside(self, r, u) result(isIt)
    class(simpleCell), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i
    logical(defBool)                        :: halfspace, sense

    ! Keep compiler happy (in immpossible case of cell with no surfaces)
    isIt = .false.

    do i= 1, size(self % surfaces)
      sense = self % surfaces(i) % surfIdx > 0
      halfspace = self % surfaces(i) % ptr % halfspace(r, u)

      isIt = halfspace .eqv. sense

      ! If halfspace is not equivalent to sense it means that point
      ! is outside the cell.
      if(.not.isIt) return
    end do

  end function inside

  !!
  !! Return distance to cell boundary
  !!
  !! See cell_inter for details
  !!
  pure subroutine distance(self, d, surfIdx, r, u)
    class(simpleCell), intent(in)           :: self
    real(defReal), intent(out)              :: d
    integer(shortInt), intent(out)          :: surfIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i
    real(defReal)                           :: test_d

    d = INF
    surfIdx = 0

    do i = 1, size(self % surfaces)
      test_d = self % surfaces(i) % ptr % distance(r, u)

      ! Select minimum distance
      if (test_d < d) then
        d = test_d
        surfIdx = i
      end if
    end do

  end subroutine distance

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(simpleCell), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Clean local
    if(allocated(self % surfaces)) deallocate (self % surfaces)

  end subroutine kill

end module simpleCell_class
