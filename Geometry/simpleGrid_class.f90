module simpleGrid_class

  use numPrecision
  use universalVariables,    only : surface_tol
  use genericProcedures,     only : fatalError, numToChar
  use dictionary_class,      only : dictionary
  use geometry_inter,        only : geometry
  use dynArray_class,        only : dynIntArray
  use nuclearDatabase_inter, only : nuclearDatabase
  use particle_class,        only : particle

  !!
  !!
  !!
  type, private :: gridCell
    integer(shortInt), dimension(:), allocatable :: mats
    real(defReal)                                :: majorant

  end type gridCell

  !!
  !! As in latUniverse_class, idx is 1 in bottom X, Y & Z corner.
  !! It increases first with X then Y and lastly Z.
  !!
  !! sizeN  -> array [nx, ny, nz], the dimensions of the grid
  !! dx     -> array [dx, dy, dz], the discretisation in each direction
  !! bounds -> [x_min, y_min, z_min, z_max, y_max, z_max] as in geometry_inter
  !!
  type, public :: grid
    class(geometry), pointer                     :: mainGeom => null()
    class(nuclearDatabase), pointer              :: xsData   => null()
    integer(shortInt), dimension(:), allocatable :: sizeN
    real(defReal), dimension(3)                  :: dx = 0
    real(defReal), dimension(6)                  :: bounds
    real(defReal), dimension(3)                  :: corner
    type(gridCell), dimension(:), allocatable    :: gridCells

  contains
    procedure :: init
    !procedure :: kill
    procedure :: getDistance
    procedure :: getValue
    procedure :: storeMats
    procedure :: update
 
  end type grid

contains

  subroutine init(self, dict, geom, xsData)
    class(grid), intent(inout)                   :: self
    class(dictionary), intent(in)                :: dict
    class(geometry), intent(in), pointer         :: geom
    class(nuclearDatabase), intent(in), pointer  :: xsData
    integer(shortInt)                            :: N
    integer(shortInt), dimension(:), allocatable :: searchN

    ! Store pointer to main geometry and data
    self % mainGeom => geom
    self % xsData   => xsData

    ! Store settings
    call dict % get(self % sizeN, 'dimensions')
    call dict % get(searchN, 'searchN')

    ! Get bounds of grid and calculate discretisations
    self % bounds = geom % bounds()

    self % dx(1)  = (self % bounds(4) - self % bounds(1)) / self % sizeN(1)
    self % dx(2)  = (self % bounds(5) - self % bounds(2)) / self % sizeN(2)
    self % dx(3)  = (self % bounds(6) - self % bounds(3)) / self % sizeN(3)

    self % corner = [self % bounds(1), self % bounds(2), self % bounds(3)]

    ! Allocate space for cells
    N = self % sizeN(1) * self % sizeN(2) * self % sizeN(3)
    allocate(self % gridCells(N))

    ! Find material idxs present in each cell
    call self % storeMats(searchN)

  end subroutine init


  !!
  !! May have issues with non-box geometry root universe surface with reflective boundary
  !!
  function getDistance(self, r, u) result(dist)
    class(grid), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: dist
    real(defReal), dimension(3)             :: rbar, point, corner, ratio
    character(100), parameter :: Here = 'getDistance (simpleGrid_class.f90)'

    ! Calculate position in grid
    rbar = r - self % corner
    point = rbar / self % dx

    ! Round each dimension either up or down depending on which boundary will be hit
    do i = 1, 3
      if (u(i) >= 0) then
        corner(i) = ceiling(point(i))
      else
        corner(i) = floor(point(i))
      end if
      ! Adjust if starting position was on boundary
      if (corner(i) == point(i)) then
        corner(i) = corner(i) + sign(ONE, u(i))
      end if
    end do

    ! Convert back to spatial coordinates - this is now the coordinates of the corner being travelled towards
    corner = corner * self % dx

    ! Determine which axis boundary will be hit first
    ratio = (corner - rbar) / u

    dist = minval(ratio)

    if (dist <= ZERO) call fatalError(Here, 'Distance invalid: '//numToChar(dist))

  end function getDistance


  !!
  !! Returns value of grid cell at position
  !!
  function getValue(self, r, u) result(val)
    class(grid), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: val
    real(defReal), dimension(3)             :: rbar
    integer(shortInt), dimension(3)         :: corner
    character(100), parameter :: Here = 'get (simpleGrid_class.f90)'

    ! Get grid cell bottom corner
    rbar = r - self % corner
    corner = floor(rbar)
    do i = 1, 3
      if (corner(i) == rbar(i) .and. u(i) < 0) then
        ! Adjust for point starting on cell boundary
        corner(i) = corner(i) - 1
      end if
    end do

    ! Adjust for bottom corner starting at 1
    corner = corner + 1

    ! Get grid cell idx
    idx = get_idx(corner, self % sizeN)
    if (idx == 0) call fatalError(Here, 'Point is outside lattice')

    val = self % gridCells(idx) % majorant

  end function getValue

  !!
  !!
  !!
  subroutine storeMats(self, searchN)
    class(grid), intent(inout)                    :: self
    integer(shortInt), dimension(3), intent(in)   :: searchN
    real(defReal), dimension(3)                   :: searchRes
    integer(shortInt)                             :: i, j, k, l, matIdx, id
    real(defReal), dimension(3)                   :: corner, r
    type(dynIntArray)                             :: mats

    ! Calculate distance between search points
    searchRes = self % dx / (searchN + 1)

    ! Loop through grid cells
    do i = 1, size(self % gridCells)

      ! Get cell lower corner
      corner = self % dx * (get_ijk(i, self % sizeN) - 1)

      ! Loop through search locations
      do j = 1, searchN(1)
        do k = 1, searchN(2)
          do l = 1, searchN(3)
            ! Find matIdx at search location
            r = corner + [j, k, l] * searchRes
            call self % mainGeom % whatIsAt(matIdx, id, r)

            ! Add to array if not already present
            if (mats % isPresent(matIdx)) then
              ! Do nothing
            else
              call mats % add(matIdx) 
            end if

          end do
        end do
      end do

      ! Store matIdx data in grid cell
      self % gridCells(i) % mats = mats % expose()

    end do

  end subroutine storeMats

  !!
  !!
  !!
  subroutine update(self)
    class(grid), intent(inout)   :: self
    integer(shortInt)            :: i, j, matIdx
    real(defReal)                :: sigmaT
    class(particle), allocatable :: p

    allocate(p)
    p % G = 1

    ! Loop through grid cells
    do i = 1, size(self % gridCells)
      ! Reset majorant
      self % gridCells(i) % majorant = ZERO

      do j = 1, size(self % gridCells(i) % mats)
        ! Get opacity of each material
        matIdx = self % gridCells(i) % mats(j)
        sigmaT = self % xsData % getTransMatXS(p, matIdx)

        ! Update majorant if required
        if (sigmaT > self % gridCells(i) % majorant) self % gridCells(i) % majorant = sigmaT

      end do

    end do

  end subroutine update



  !!
  !! Generate ijk from localID and shape
  !!
  !! Args:
  !!   localID [in] -> Local id of the cell between 1 and product(sizeN)
  !!   sizeN [in]   -> Number of cells in each cardinal direction x, y & z
  !!
  !! Result:
  !!   Array ijk which has integer position in each cardinal direction
  !!
  pure function get_ijk(localID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: localID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = localID - 1

    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = temp - sizeN(2) * base + 1

    ijk(3) = base + 1

  end function get_ijk


  pure function get_idx(ijk, sizeN) result(idx)
    integer(shortInt), dimension(3), intent(in) :: ijk
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt)                           :: idx

    if (any(ijk <= 0 .or. ijk > sizeN)) then ! Point is outside lattice
      idx = 0
    else
      idx = ijk(1) + sizeN(1) * (ijk(2)-1 + sizeN(2) * (ijk(3)-1))
    end if

  end function get_idx

end module simpleGrid_class
