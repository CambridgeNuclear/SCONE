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
  !!    #origin (1.0 0.0);# // Optional. Default (0 0)
  !!    rGrid lin;
  !!    #Rmin 2.0;#  // Optional. Default 0.0
  !!    Rmax  10.0;
  !!    rN 10;
  !!    #zGrid lin; // Optional
  !!    Zmin 2.0;
  !!    Zmax  10.0;
  !!    zN 10;
  !!    }
  !!
  !!
  type, public, extends (tallyMap) :: cylindricalMap
    private
    real(defReal), dimension(2) :: origin = ZERO
    type(grid)                  :: rBounds
    type(grid)                  :: zBounds
    integer(shortInt)           :: rN = 0
    integer(shortInt)           :: zN = 0
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
    real(defReal)                            :: Rmin, Rmax, Zmin, Zmax, vol
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'init (sphericalmap_class.f90)'

    ! Check & load origin
    call dict % getOrDefault(temp, 'origin', [ZERO, ZERO])

    if (size(temp) /= 2) then
      call fatalError(Here, 'Expected 2 values for origin. Got: ' // numToChar(size(temp)))
    end if
    self % origin = temp

    ! Load radial grid information
    if (.not.dict % isPresent('rGrid')) call fatalError(Here, 'Keyword rGrid must be present')
    call dict % get(type, 'rGrid')

    select case(type)
      case('lin')
        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % rN, 'N')

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
        call dict % get(self % rN, 'N')

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
    if (dict % isPresent('zGrid')) then

      call dict % get(type, 'zGrid')

      select case(type)
        case('lin')
          call dict % get(Zmin, 'Zmin')
          call dict % get(Zmax, 'Zmax')
          call dict % get(self % zN, 'N')

          ! Build grid
          call self % zBounds % init(Zmin, Zmax, self % zN, type)

        case default
          call fatalError(Here, "'zGrid' can take only values of: lin")

        end select

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
      if (self % zN /= 0) N = N * self % zN
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
    class(cylindricalMap), intent(in) :: self
    character(nameLen)          :: name

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
    integer(shortInt)                 :: idx, zIdx
    real(defReal)                     :: r

    ! Calculate the distance from the origin
    r = norm2(state % r(1:2) - self % origin)

    ! Search and return 0 if r is out-of-bounds
    idx = self % rBounds % search(r)

    if (idx == valueOutsideArray) then
      idx = 0
      return
    end if

    ! Return if there is no axial dimension
    if (self % zN == 0) return

    ! Search along the axial dimension
    zIdx = self % zBounds % search(state % r(3))
    if (zIdx == valueOutsideArray) then
      idx = 0
      return
    end if

    idx = self % rN * (zIdx - 1) + idx

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
    name = trim(self % getAxisName()) //'radialBounds'

    call out % startArray(name, [2,self % rN])
    do i = 1, self % rN
      ! Print lower bin boundary
      call out % addValue(self % rBounds % bin(i))

      ! Print upper bin boundar
      call out % addValue(self % rBounds % bin(i + 1))

    end do
    call out % endArray()

    ! Add axial dimension if present
    if (self % zN /= 0) then
      name = trim(self % getAxisName()) //'axialBounds'

      call out % startArray(name, [2,self % zN])
      do i = 1, self % zN
        ! Print lower bin boundary
        call out % addValue(self % zBounds % bin(i))

        ! Print upper bin boundar
        call out % addValue(self % zBounds % bin(i + 1))

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
    call self % zBounds % kill()
    self % rN = 0
    self % zN = 0

  end subroutine kill

end module cylindricalMap_class
