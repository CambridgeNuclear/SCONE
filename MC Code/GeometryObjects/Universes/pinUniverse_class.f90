module pinUniverse_class

  use numPrecision
  use universalVariables
  use genericProcedures,   only : fatalError, hasDuplicates, linFind, targetNotFound, swap
  use vector_class,        only : vector
  use dictionary_class,    only : dictionary
  use maps_class,          only : intMap

  use coord_class,         only : coord
  use surface_inter,       only : surface, surfaceSlot, surfaceShelf
  use surfaceFactory_func, only : new_surface
  use cell_class,          only : cell, cellShelf
  use universe_inter,      only : universe

  !*** STAYS HERE ONLY PROVISIONALLY
  use nuclearData_inter, only : nuclearData

  implicit none
  private


  !!
  !! Small private type to group all data related to given annulus
  !!
  type, private :: annulus
    integer(shortInt) :: cellIdx
    integer(shortInt) :: surfIdx
  end type annulus

  !!
  !!
  !!
  type, public :: pinUniverse
    private
    real(defReal), dimension(:), allocatable :: r
    type(annulus), dimension(:), allocatable :: data


  contains
    procedure :: init

    procedure, private :: buildAnnuli
    procedure, private :: makeFillVector
  end type pinUniverse

contains

  !!
  !! Initialise pinUniverse
  !! use u<####> syntax for the universe fill
  !! radius = 0.0 specifies outermost annulus !
  !!
  subroutine init(self, fillVector, id, matNames, radii, sShelf, cShelf, materials)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(pinUniverse), intent(inout)                      :: self
    integer(shortInt), intent(in)                          :: id
    character(nameLen),dimension(:),intent(in)             :: matNames
    real(defReal), dimension(:), intent(in)                :: radii
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(cellShelf),intent(inout)                          :: cShelf
    class(nuclearData), intent(in)                         :: materials
    character(nameLen),dimension(size(radii))              :: names
    integer(shortInt)                                      :: N, i, idx, surfID
    type(dictionary)                                       :: tempDict
    class(surface), allocatable                            :: surf
    character(100), parameter :: Here = 'init (pinUniverse_class.f90)'

    ! Check that radii are -ve and that matNames and radii are the same size
    if (any (radii < ZERO)) call fatalError(Here,'-ve radii cannot be provided')
    if (size(radii) /= size(matNames)) call fatalError(Here,'Size of matNames & radii does not match')
    if (hasDuplicates(radii)) call fatalError(Here,'Duplicate values in annuli radii')

    ! Load number of annuli
    N = size(radii)

    ! Allocate storage space and set offset to 0
    !call self % setOffset([ZERO, ZERO, ZERO])
    self % r = radii
    allocate(self % data(N))

    ! Copy material names to local copy
    names = matNames

    ! Radii may not be in right order. Sort them together with corresponding material names

    ! Start by processing the outermost element
    idx = linFind(self % r, ZERO)
    if(idx == targetNotFound) call fatalError(Here,'Outermost element with radius 0.0 was not found')

    call swap( self % r(idx),    self % r(N))
    call swap( names(idx), names(N))

    ! Put the rest in the ascending order -> Selection sort for simplicity.
    ! This is not performace critical (not yet at least)
    do i = N-1,1,-1
      idx = maxloc( self % r(1:i), 1 )
      call swap (self % r(idx), self % r(i))
      call swap (names(idx), names(i) )
    end do

    ! Build surfaces and cells
    call self % buildAnnuli(sShelf, cShelf)

    ! Create build vector
    call self % makeFillVector(fillVector, matNames, materials)


  end subroutine init

  !!
  !! Build all surfaces and cells that define pinUniverse
  !!
  subroutine buildAnnuli(self, sShelf, cShelf)
    class(pinUniverse), intent(inout)           :: self
    type(surfaceShelf), intent(inout)           :: sShelf
    type(cellShelf), intent(inout)              :: cShelf
    type(dictionary)                            :: tempDict
    class(surface),allocatable                  :: tempSurf
    type(cell)                                  :: tempCell
    integer(shortInt)                           :: cellId, i
    integer(shortInt),dimension(size(self % r)) :: surfIDs

    ! Build all surfaces
    do i =1,size(self % r)
      ! Initialise temporary dictionary
      call tempDict % init(4,1)

      ! Build approperiate surface
      if (self % r(i) == ZERO ) then
        call tempDict % store('type','infSurf')
        call tempDict % store('id', 1 )

      else
        call tempDict % store('type','zCylinder')
        call tempDict % store('id',1)
        call tempDict % store('origin',[ZERO, ZERO, ZERO])
        call tempDict % store('radius', self % r(i))

      end if

      ! Put surface on the shelf
      allocate(tempSurf, source = new_surface(tempDict))
      call sShelf % getOrAdd(tempSurf, surfIDs(i), self % data(i) % surfIdx)
      call tempDict % kill()
    end do

    ! Build innermost cell
    call tempCell % init( [-surfIDs(1)], 1, sShelf)
    call cShelf % getOrAdd(tempCell, cellId, self % data(1) % cellIdx)

    ! Build rest of the cells
    do i=2,size(self % r)
      call tempCell % init([surfIDs(i-1), -surfIDs(i)], 1, sShelf)
      call cShelf % getOrAdd(tempCell, cellId, self % data(i) % cellIdx)

    end do
  end subroutine buildAnnuli

  !!
  !! Translate material names to fillVector
  !!
  subroutine makeFillVector(self, fillVEctor, matNames ,materials)
    class(pinUniverse), intent(inout)                      :: self
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    character(nameLen),dimension(:),intent(in)             :: matNames
    class(nuclearData), intent(in)                         :: materials
    integer(shortInt)                                      :: i, fill, L
    character(nameLen)                                     :: tempName
    logical(defBool)                                       :: isUniverse
    character(100), parameter :: Here ='makeFillVector (pinUniverse_class.f90)'

    ! Double check that size of radii match size of matNames
    if( size(self % r) /= size(matNames)) call fatalError(Here,'matNames not matching self % r')

    ! Allocate fill vector
    allocate( fillVector( size(matNames)))

    ! Loop over all matNames
    do i=1, size(matNames)
      ! Copy material name, adjust left and find length of name
      tempName = adjustl(matNames(i))
      L = len_trim(tempName)

      ! Verify if name is a universe idd
      if ( tempName(1:2) == 'u<' .and. tempName(L:L) == '>') then
        ! Read universe id into fill
        tempName = tempName(3:L-1)
        read(tempName,*) fill
        fill = -fill

      else
        fill = materials % getIdx(tempName)

      end if

       fillVector(i) = fill
    end do

      print *, fillVector

  end subroutine makeFillVector


end module pinUniverse_class
