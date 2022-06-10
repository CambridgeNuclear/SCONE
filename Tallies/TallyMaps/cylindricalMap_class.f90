module cylindricalMap_class

  use numPrecision
  use universalVariables, only : valueOutsideArray
  use genericProcedures,  only : fatalError, dotProduct, numToChar
  use dictionary_class,   only : dictionary
  use grid_class,         only : grid
  use particle_class,     only : particleState
  use outputFile_class,   only : outputFile
  use tallyMap_inter,     only : tallyMap, kill_super => kill

  implicit none
  private


  !!
  !! Divides space into a mesh in cylindrical co-ordinates
  !!
  !! TODO:
  !!   Implement polar & azimuthal subdivision
  !!   Add unstructured grid boundaries
  !!
  !! Interface:
  !!   tallyMap interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  cylindricalMap {
  !!    type cylindricalMap;
  !!    #orientation x;# // Optional. Defeult z
  !!    #origin (1.0 0.0);# // Optional. Default (0 0). Order is xy, xz or yz.
  !!    rGrid lin;
  !!    #Rmin 2.0;#  // Optional. Default 0.0
  !!    Rmax  10.0;
  !!    rN 10;
  !!    #axGrid lin; // Optional
  !!    #axMin 2.0;   // Needed if axGrid is present
  !!    #axMax  10.0; // Needed if axGrid is present
  !!    #axN 10;      // Needed if axGrid is present
  !!    #azimuthalBins 4;# // Optional
  !!    }
  !!
  !!
  type, public, extends (tallyMap) :: cylindricalMap
    private
    real(defReal), dimension(2) :: origin = ZERO
    type(grid)                  :: rBounds
    type(grid)                  :: axBounds
    type(grid)                  :: azBounds
    integer(shortInt)           :: rN = 0
    integer(shortInt)           :: axN = 0
    integer(shortInt)           :: azN = 0
    integer(shortInt)           :: DIM1
    integer(shortInt)           :: DIM2
    integer(shortInt)           :: DIM3
  contains
    ! Superclass
    procedure :: init
    procedure :: bins
    procedure :: dimensions
    procedure :: getAxisName
    procedure :: map
    procedure :: print
    procedure :: kill
  end type cylindricalMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(cylindricalMap), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: temp, grid
    character(nameLen)                       :: type
    real(defReal)                            :: Rmin, Rmax, Amin, Amax, vol
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'init (sphericalmap_class.f90)'

    ! Check & load origin
    call dict % getOrDefault(temp, 'origin', [ZERO, ZERO])

    if (size(temp) /= 2) then
      call fatalError(Here, 'Expected 2 values for origin. Got: ' // numToChar(size(temp)))
    end if
    self % origin = temp

    ! Check orientation of the cylinder
    if (dict % isPresent('orientation')) then
      select case(type)
        case('x')
          self % DIM1 = 2
          self % DIM2 = 3
          self % DIM3 = 1
        case('y')
          self % DIM1 = 1
          self % DIM2 = 3
          self % DIM3 = 2
        case('z')
          self % DIM1 = 1
          self % DIM2 = 2
          self % DIM3 = 3
      end select
    else
      self % DIM1 = 1
      self % DIM2 = 2
      self % DIM3 = 3
    end if

    ! Load radial grid information
    if (.not.dict % isPresent('rGrid')) call fatalError(Here, 'Keyword rGrid must be present')
    call dict % get(type, 'rGrid')

    ! Check type of radial grid bins
    select case(type)
      case('lin')
        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % rN, 'rN')

        ! Check that minimum radius is OK
        if (Rmin < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve is:'//numToChar(Rmin))
        end if

        ! Build grid
        call self % rBounds % init(Rmin, Rmax, self % rN, type)

      case('unstruct')

        call dict % get(grid,'bins')

        ! Initialise
        self % rN = size(grid) - 1
        call self % rBounds % init(grid)

      case('equivolume')

        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % rN, 'rN')

        ! Check that minimum radius is OK
        if (Rmin < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve is:'//numToChar(Rmin))
        end if

        ! Calculate volume
        vol = (Rmax**2 - Rmin**2) /self % rN

        allocate(grid(self % rN + 1))
        ! Calculate grid boundaries
        grid(1) = Rmin
        do i = 2,size(grid)
          grid(i) = sqrt(vol + grid(i-1)**2)
        end do

        call self % rBounds % init(grid)

      case default
        call fatalError(Here, "'rGrid' can take only values of: lin")

    end select

    ! Load axial grid information
    if (dict % isPresent('axGrid')) then

      call dict % get(type, 'axGrid')

      select case(type)
        case('lin')
          call dict % get(Amin, 'axMin')
          call dict % get(Amax, 'axMax')
          call dict % get(self % axN, 'axN')

          ! Build grid
          call self % axBounds % init(Amin, Amax, self % axN, type)

        case default
          call fatalError(Here, "'axGrid' can take only values of: lin")

        end select

      end if

      ! Load azimuthal grid information
      if (dict % isPresent('azimuthalBins')) then
        call dict % get(self % azN, 'azimuthalBins')
        call self % azBounds % init(-PI, PI, self % azN, 'lin')
      end if

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(cylindricalMap), intent(in) :: self
    integer(shortInt), intent(in)     :: D
    integer(shortInt)                 :: N

    if (D == 1 .or. D == 0) then
      N = self % rN
      if (self % axN /= 0) N = N * self % axN
      if (self % azN /= 0) N = N * self % azN
    else
      N = 0
    end if

  end function bins

  !!
  !! Return number of dimensions
  !!
  !! See tallyMap for specification.
  !!
  elemental function dimensions(self) result(D)
    class(cylindricalMap), intent(in)  :: self
    integer(shortInt)                  :: D

    D = 1

  end function dimensions

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(cylindricalMap), intent(in)  :: self
    character(nameLen)                 :: name

    name ='cylindricalMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self, state) result(idx)
    class(cylindricalMap), intent(in) :: self
    class(particleState), intent(in)  :: state
    integer(shortInt)                 :: idx, aIdx, N
    real(defReal)                     :: r, x, y, theta

    ! Calculate the distance from the origin
    r = norm2(state % r([self % DIM1, self % DIM2]) - self % origin)

    ! Search and return 0 if r is out-of-bounds
    idx = self % rBounds % search(r)
    N = self % rN

    if (idx == valueOutsideArray) then
      idx = 0
      return
    end if

    ! Update if there is axial dimension
    if (self % axN /= 0) then
      ! Search along the axial dimension and return 0 if index is out-of-bounds
      aIdx = self % axBounds % search(state % r(self % DIM3))
      if (aIdx == valueOutsideArray) then
        idx = 0
        return
      end if
      ! Compute new index and update number of bins
      idx = N * (aIdx - 1) + idx
      N = N * self % axN
    end if

    ! Update if there is azimuthal dimension
    if (self % azN /= 0) then

      x = state % r(self % DIM1) - self % origin(1)
      y = state % r(self % DIM2) - self % origin(2)

      theta = atan2(y,x)
      ! Search along the azimuthal dimension and return 0 if index is out-of-bounds
      aIdx = self % azBounds % search(theta)
      ! Compute new index and update number of bins
      idx = N * (aIdx - 1) + idx
    end if

  end function map

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification.
  !!
  subroutine print(self,out)
    class(cylindricalMap), intent(in)  :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'RadialBounds'

    call out % startArray(name, [2,self % rN])
    do i = 1, self % rN
      ! Print lower bin boundary
      call out % addValue(self % rBounds % bin(i))

      ! Print upper bin boundar
      call out % addValue(self % rBounds % bin(i + 1))

    end do
    call out % endArray()

    ! Add axial dimension if present
    if (self % axN /= 0) then
      name = trim(self % getAxisName()) //'AxialBounds'

      call out % startArray(name, [2,self % axN])
      do i = 1, self % axN
        ! Print lower bin boundary
        call out % addValue(self % axBounds % bin(i))

        ! Print upper bin boundar
        call out % addValue(self % axBounds % bin(i + 1))

      end do
      call out % endArray()
    end if

    ! Add azimuthal dimension if present
    if (self % azN /= 0) then
      name = trim(self % getAxisName()) //'AzimuthalBounds'

      call out % startArray(name, [2,self % azN])
      do i = 1, self % azN
        ! Print lower bin boundary
        call out % addValue(self % azBounds % bin(i) + PI)

        ! Print upper bin boundar
        call out % addValue(self % azBounds % bin(i + 1) + PI)

      end do
      call out % endArray()
    end if

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cylindricalMap), intent(inout) :: self

    call kill_super(self)

    self % origin = ZERO
    call self % rBounds % kill()
    call self % axBounds % kill()
    call self % azBounds % kill()
    self % rN = 0
    self % axN = 0
    self % azN = 0

  end subroutine kill

end module cylindricalMap_class
