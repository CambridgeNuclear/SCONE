!
! Uniformly spaced 2D/3D rectangular mesh
!
module mesh_class
  use numPrecision
  use genericProcedures

  implicit none

  type, public :: mesh
    real(defReal), dimension(3) :: corner                        ! bottom corner of the mesh
    integer(shortInt), dimension(3) :: extent                    ! the number of mesh cells in each direction
    real(defReal), dimension(3) :: pitch                         ! width of mesh cells in each direction
    logical(defBool) :: is3D                                     ! is the mesh 3D?
    integer(longInt),dimension(:,:,:), allocatable :: density    ! stores the density of, say, collision (may be moved to tallies later!)
    integer(shortInt) :: id                                      ! unique mesh id
    character(100), public :: name = ""
  contains
    procedure :: init             ! initialise mesh
    procedure :: whichMeshCell    ! Identify which cell a point occupies
    procedure :: score            ! Score an additional 'hit' in a mesh cell
  end type mesh

contains

  subroutine init(self, pitch, corner, sz, id, is3D, name)
    class(mesh), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: pitch
    real(defReal), dimension(3), intent(in) :: corner
    integer(shortInt), dimension(3), intent(in) :: sz
    integer(shortInt), intent(in) :: id
    character(*), intent(in), optional :: name
    logical(defBool), intent(in) :: is3D

    self % pitch = pitch
    self % extent = sz
    self % corner = corner
    self % is3D = is3D
    self % id = id
    if(is3D) then
      allocate(self % density(sz(1),sz(2),sz(3)))
    else
      allocate(self % density(sz(1),sz(2),1))
      self % density(:,:,:) = 0_longInt
    end if
    if(present(name)) self % name = name
  end subroutine init

  !
  ! Return which cell in the mesh a point occupies
  ! Indexing goes from bottom left corner to top right
  !
  function whichMeshCell(self,r,u)result(ijk)
    class(mesh), intent(in) :: self
    real(defReal), intent(in), dimension(3) :: r
    real(defReal), intent(in), dimension(3) :: u
    real(defReal), dimension(3) :: corner, pitch
    integer(shortInt), dimension(3) :: ijk

    corner = self % corner
    pitch = self % pitch

    ijk(1) = ceiling((r(1) - corner(1))/pitch(1))
    ijk(2) = ceiling((r(2) - corner(2))/pitch(2))
    if (self % is3D) then
      ijk(3) = ceiling((r(3) - corner(3))/pitch(3))
    else
      ijk(3) = 1
    end if
  end function whichMeshCell

  subroutine score(self,r,u)
    class(mesh), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r,u
    integer(shortInt), dimension(3) :: ijk
    logical(defBool) :: inMesh

    ijk = self % whichMeshCell(r,u)
    inMesh = all((ijk>0).and.(ijk<=self%extent))
    if(inMesh) self % density(ijk(1),ijk(2),ijk(3)) = self % density(ijk(1),ijk(2),ijk(3)) + 1
  end subroutine score

end module mesh_class
