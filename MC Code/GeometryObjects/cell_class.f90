!
! The cell class is the basic geometrical structure: it contains either universes/lattices or materials
! This class will contains bounding surfaces and their halfspace and can be searched to identify
! whether a particle is inside the cell.
!
module cell_class
  use numPrecision
  use genericProcedures
  use universalVariables
  use surface_class

  implicit none
  private

  type, public :: cell
    type(surface_ptr), dimension(:), allocatable :: surfaces ! the surfaces which define the cell
    logical(defBool), dimension(:), allocatable :: halfspaces ! The halfspace of each surface corresponding to inside the cell
    integer(shortInt) :: numSurfaces                          ! the number of surfaces which define the cell
    integer(shortInt) :: fillType = 0                         ! determines if cell contains a material, universe, or lattice (1,2,3)
    integer(shortInt) :: uniInd = 0                           ! index of the universe filling the cell
    integer(shortInt) :: latInd = 0                           ! index of the lattice filling the cell
    integer(shortInt) :: materialInd = 0                      ! index of the material contained in cel
    integer(shortInt) :: id                                   ! unique user-defined ID
    integer(shortInt) :: parentUni = 0                        ! index of the parent universe of the cell
    logical(defBool) :: insideGeom = .true.                   ! is cell within geometry? Used to invoke BCs
    integer(shortInt) :: instances = 1                        ! the number of instances of a given cell
    integer(shortInt) :: geometryInd = 0                      ! the index of the cell in the full geometry array
    real(defReal) :: volume                                   ! the volume of the cell
    character(100), public :: name = ""
  contains
    procedure :: init
    procedure :: fill
    procedure :: insideCell
    procedure :: getDistance
    procedure :: whichSurface
  end type cell

  ! Wrapper type for cell pointers
  type, public :: cell_ptr
    class(cell), pointer :: ptr
  contains
    procedure :: init => init_ptr
    procedure :: fill => fill_ptr
    procedure :: insideCell => insideCell_ptr
    procedure :: getDistance => getDistance_ptr
    procedure :: whichSurface => whichSurface_ptr
    procedure :: insideGeom => insideGeom_ptr
    procedure :: fillType => fillType_ptr
    procedure :: uniInd => uniInd_ptr
    procedure :: latInd => latInd_ptr
    procedure :: materialInd => materialInd_ptr
    procedure :: geometryInd => geometryInd_ptr
    procedure :: associated => associated_ptr
    procedure :: name => name_ptr
    procedure :: kill
    generic   :: assignment(=) => cell_ptr_assignment, cell_ptr_assignment_target
    procedure,private :: cell_ptr_assignment
    procedure,private :: cell_ptr_assignment_target
    !procedure,private :: cell_ptr_assignment_null
  end type cell_ptr

contains

  !!
  !! Initialise the cell with the bounding surfaces, their halfspaces and parent universe
  !!
  subroutine init(self, surfaces, halfspaces, id, fillType, geometryInd, name)
    class(cell), intent(inout) :: self
    class(surface_ptr), dimension(:), intent(in) :: surfaces
    logical(defBool), dimension(:), intent(in) :: halfspaces
    integer(shortInt), intent(in) :: id
    integer(shortInt), intent(in) :: fillType
    integer(shortInt), intent(in) :: geometryInd
    character(*), optional, intent(in) :: name

    self % numSurfaces = size(halfspaces)
    allocate(self % surfaces(self % numSurfaces))
    allocate(self % halfspaces(self % numSurfaces))
    self % surfaces = surfaces
    self % halfspaces = halfspaces
    self % id = id
    self % fillType = fillType
    self % geometryInd = geometryInd
    ! If outside of the geometry, parentUni will be meaningless
    if (fillType == outsideFill) then
      self % insideGeom = .false.
    end if
    if (present(name)) self % name = name

  end subroutine init

  !!
  !! Fills the cell with a material, universe or lattice
  !!
  subroutine fill(self, fillInd)
    class(cell), intent(inout) :: self
    integer(shortInt), intent(in) :: fillInd

    ! Fill cell depending on given fill type
    if (self % fillType == materialFill) then
      self % materialInd = fillInd
    else if (self % fillType == universeFill) then
      self % uniInd = fillInd
    else if (self % fillType == latticeFill) then
      self % latInd = fillInd
    else if (self % fillType /= outsideFill) then
      call fatalError('fill, cell', 'Cell filled with an incorrect fill index')
    end if
  end subroutine fill

  !!
  !! Checks whether a point occupies a cell by examining each surface in turn
  !! Returns TRUE if point is in cell, FALSE if point is not
  !!
  function insideCell(self, r, u) result(exists)
    class(cell), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    logical(defBool) :: exists, &  ! whether point exists in cell
                        sense      ! halfspace of surface in which point exists
    integer(shortInt) :: i

    ! Need only check that halfspaces are satisfied: if not, the point is outside the cell
    ! Check each surrounding surface of the cell to ensure point is within its halfspace
    do i = 1, self % numSurfaces
      exists = self % surfaces(i) % halfspace(r,u)
      sense = self % halfspaces(i)

      ! If not in halfspace, terminate the search
      if ( exists .neqv. sense ) then
        exists = outside
        return
      end if
    end do
    exists = inside
    return

  end function insideCell

  !!
  !! Find the shortest positive distance to the boundary of the cell
  !!
  function getDistance(self, r, u) result(distance)
    class(cell), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal) :: distance
    real(defReal) :: testDistance
    integer(shortInt) :: i

    distance = INFINITY

    ! Search through all surfaces to find the minimum distance to the boundary
    ! Should not have to check for negative distances: these are set to INFINITY
    ! in the distance routines
    do i = 1, self%numSurfaces
      testDistance = self % surfaces(i) % distanceToSurface(r,u)
      if (testDistance < distance) then
        distance = testDistance
      end if
    end do

  end function getDistance

  !!
  !! Find which surface of a cell was crossed by a particle
  !! For compound surfaces must return component surface
  !!
  function whichSurface(self, r, u)result(surfPointer)
    class(cell), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    integer(shortInt) :: i
    real(defReal) :: distance, testDistance

    distance = INFINITY
    ! First identify which surface will have been crossed assuming all surfaces are simple
    do i = 1, self % numSurfaces
      testDistance = self % surfaces(i) % distanceToSurface(r,u)
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % surfaces(i) % ptr
      end if
    end do

    ! If the identified surface is compound, identify which constituent surface is crossed
    if (surfPointer % isCompound) then
      surfPointer => surfPointer % whichSurface(r, u)
    end if


  end function whichSurface

!
! Pointer wrapper functions
!
  !!
  !! Initialise the cell with the bounding surfaces, their halfspaces and the
  !! indices describing the fill type, fill, and parent universe
  !!
  subroutine init_ptr(self, surfaces, halfspaces, id, fillType, geometryInd, name)
    class(cell_ptr), intent(inout) :: self
    class(surface_ptr), dimension(:), intent(in) :: surfaces
    logical(defBool), dimension(:), intent(in) :: halfspaces
    integer(shortInt), intent(in) :: id
    integer(shortInt), intent(in) :: fillType
    integer(shortInt), intent(in) :: geometryInd
    character(100), optional, intent(in) :: name
    call self % ptr % init(surfaces, halfspaces, id, fillType, geometryInd, name)
  end subroutine init_ptr

  subroutine fill_ptr(self, fillInd)
    class(cell_ptr), intent(inout) :: self
    integer(shortInt), intent(in) :: fillInd
    call self % ptr % fill(fillInd)
  end subroutine fill_ptr

  !!
  !! Checks whether a point occupies a cell by examining each surface in turn
  !! Returns TRUE if point is in cell, FALSE if point is not
  !!
  function insideCell_ptr(self, r, u) result(exists)
    class(cell_ptr), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool) :: exists
    exists = self % ptr % insideCell(r,u)
  end function insideCell_ptr

  !!
  !! Find the shortest positive distance to the boundary of the cell
  !!
  function getDistance_ptr(self, r, u) result(distance)
    class(cell_ptr), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal) :: distance
    distance = self % ptr % getDistance(r,u)
  end function getDistance_ptr

  !!
  !! Return whether cell_ptr points to a cell which is inside the geometry
  !!
  function insideGeom_ptr(self)result(insideGeom)
    class(cell_ptr), intent(in) :: self
    logical(defBool) :: insideGeom
    insideGeom = self % ptr % insideGeom
  end function insideGeom_ptr

  !!
  !! Returns the fill type of the cell pointed to by cell_ptr
  !!
  function fillType_ptr(self)result(fillType)
    class(cell_ptr), intent(in) :: self
    integer(shortInt) :: fillType
    fillType = self % ptr % fillType
  end function fillType_ptr

  !!
  !! Returns the universe index of the cell pointed to by cell_ptr
  !!
  function uniInd_ptr(self)result(uniInd)
    class(cell_ptr), intent(in) :: self
    integer(shortInt) :: uniInd
    uniInd = self % ptr % uniInd
  end function uniInd_ptr

  !!
  !! Returns the lattice index of the cell pointed to by cell_ptr
  !!
  function latInd_ptr(self)result(latInd)
    class(cell_ptr), intent(in) :: self
    integer(shortInt) :: latInd
    latInd = self % ptr % latInd
  end function latInd_ptr

  !!
  !! Returns the material index of the cell pointed to by cell_ptr
  !!
  function materialInd_ptr(self)result(materialInd)
    class(cell_ptr), intent(in) :: self
    integer(shortInt) :: materialInd
    materialInd = self % ptr % materialInd
  end function materialInd_ptr

  !!
  !! Returns the geometry index of the cell pointed to by cell_ptr
  !!
  function geometryInd_ptr(self)result(geometryInd)
    class(cell_ptr), intent(in) :: self
    integer(shortInt) :: geometryInd
    geometryInd = self % ptr % geometryInd
  end function geometryInd_ptr

  !!
  !! Check whether the pointer wrapper is associated to a cell
  !!
  function associated_ptr(self)result(assoc)
    class(cell_ptr), intent(in) :: self
    logical(defBool) :: assoc
    assoc = associated(self % ptr)
  end function associated_ptr

  !!
  !! Returns the name of the cell pointed to by cell_ptr
  !!
  function name_ptr(self)result(name)
    class(cell_ptr), intent(in) :: self
    character(100) :: name
    name = self % ptr % name
  end function name_ptr

  !!
  !! Find which surface of a cell was crossed by a particle
  !! For compound surfaces must return component surface
  !!
  function whichSurface_ptr(self, r, u)result(surfPointer)
    class(cell_ptr), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    surfPointer => self % ptr % whichSurface(r,u)
  end function whichSurface_ptr

  subroutine cell_ptr_assignment(LHS,RHS)
    class(cell_ptr), intent(out)  :: LHS
    type(cell_ptr), intent(in)    :: RHS

    !if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr
  end subroutine cell_ptr_assignment

  subroutine cell_ptr_assignment_target(LHS,RHS)
    class(cell_ptr), intent(out)        :: LHS
    class(cell), target, intent(in)     :: RHS

    !if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS

  end subroutine cell_ptr_assignment_target

  subroutine kill(self)
    class(cell_ptr), intent(inout) :: self
    self % ptr => null()
  end subroutine kill

end module cell_class
