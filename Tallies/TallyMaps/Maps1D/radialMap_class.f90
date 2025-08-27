module radialMap_class

  use numPrecision
  use universalVariables, only : valueOutsideArray, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use grid_class,         only : grid
  use particle_class,     only : particleState
  use outputFile_class,   only : outputFile
  use tallyMap1D_inter,   only : tallyMap1D, kill_super => kill

  implicit none
  private

  !!
  !! Divides space into a radial mesh in cylindrical or spherical co-ordinates
  !!
  !! Interface:
  !!   tallyMap1D interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  radialMap {
  !!    type radialMap;
  !!    axis x;                 // Optional. Default is 'xyz' (spherical map)
  !!    #origin (1.0 0.0 0.0);# // Optional. Default (0.0 0.0 0.0)
  !!    grid lin;
  !!    #min 2.0;#              // Optional. Default 0.0
  !!    max  10.0;
  !!    N 10;
  !!    }
  !!
  type, public, extends (tallyMap1D) :: radialMap
    private
    integer(shortInt), dimension(:), allocatable :: axis
    real(defReal), dimension(3) :: origin
    type(grid)                  :: bounds
    integer(shortInt)           :: N = 0

  contains
    ! Superclass
    procedure :: init
    procedure :: bins
    procedure :: getAxisName
    procedure :: map
    procedure :: print
    procedure :: kill

  end type radialMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(radialMap), intent(inout)          :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: temp, grid
    character(nameLen)                       :: type
    real(defReal)                            :: min, max, vol, exp
    integer(shortInt)                        :: i
    logical(defBool)                         :: spherical
    character(100), parameter :: Here = 'init (radialMap_class.f90)'

    ! Check if the map is cylindrical or spherical, and orientation of the cylinder
    call dict % getOrDefault(type, 'axis', 'xyz')
    spherical = .false.

    select case(type)
      case('x')
        allocate(self % axis(2))
        self % axis = [Y_AXIS, Z_AXIS]

      case('y')
        allocate(self % axis(2))
        self % axis = [X_AXIS, Z_AXIS]

      case('z')
        allocate(self % axis(2))
        self % axis = [X_AXIS, Y_AXIS]

      case('xyz')
        allocate(self % axis(3))
        self % axis = [X_AXIS, Y_AXIS, Z_AXIS]
        spherical = .true.

      case default
        call fatalError(Here, 'Keyword orientation must be x, y, z or xyz for spherical. It is: '//type)

    end select

    ! Check & load origin
    call dict % getOrDefault(temp, 'origin', [ZERO, ZERO, ZERO])

    if (size(temp) /= 3) then
      call fatalError(Here, 'Expected 3 values for origin. Got: ' // numToChar(size(temp)))
    end if

    self % origin = temp

    ! Load radial grid information
    if (.not. dict % isPresent('grid')) call fatalError(Here, 'Keyword grid must be present')
    call dict % get(type, 'grid')

    ! Check type of radial grid bins
    select case(type)

      case('lin')

        call dict % getOrDefault(min, 'min', ZERO)
        call dict % get(max, 'max')
        call dict % get(self % N, 'N')

        ! Check that minimum radius is ok
        if (min < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve. It is: '//numToChar(min))
        end if

        ! Build grid
        call self % bounds % init(min, max, self % N, type)

      case('unstruct')

        call dict % get(grid,'bins')

        ! Initialise
        self % N = size(grid) - 1
        call self % bounds % init(grid)

      case('equivolume')

        call dict % getOrDefault(min, 'min', ZERO)
        call dict % get(max, 'max')
        call dict % get(self % N, 'N')

        ! Check that minimum radius is OK
        if (min < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve. It is: '//numToChar(min))
        end if

        ! Get exponent  for spherical or cylindrical mesh
        if (spherical) then
          exp = 3.0_defReal
        else
          exp = 2.0_defReal
        end if

        ! Allocate grid and initialise grid boundaries
        allocate(grid(self % N + 1))
        grid(1) = min
        grid(self % N + 1) = max

        ! Calculate volume
        vol = (max**exp - min**exp)/self % N

        ! Calculate grid boundaries
        do i = 2,self % N
          grid(i) = (vol + grid(i-1)**exp)**(ONE/exp)
        end do

        call self % bounds % init(grid)

      case default
        call fatalError(Here, "'grid' can take only values of: lin, unstruct, equivolume")

    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(radialMap), intent(in) :: self
    integer(shortInt), intent(in)        :: D
    integer(shortInt)                    :: N

    if (D == 1 .or. D == 0) then
      N = self % N
    else
      N = 0
    end if

  end function bins

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(radialMap), intent(in) :: self
    character(nameLen)                   :: name

    name = 'radialMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self, state) result(idx)
    class(radialMap), intent(in)     :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx
    real(defReal)                    :: r

    ! Calculate the distance from the origin
    r = norm2(state % r(self % axis) - self % origin(self % axis))

    ! Search and return 0 if r is out-of-bounds
    idx = self % bounds % search(r)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification.
  !!
  subroutine print(self,out)
    class(radialMap), intent(in) :: self
    class(outputFile), intent(inout)     :: out
    character(nameLen)                   :: name
    integer(shortInt)                    :: i

    ! Name the array
    name = trim(self % getAxisName()) //'RadialBounds'

    call out % startArray(name, [2,self % N])
    do i = 1, self % N
      ! Print lower bin boundary
      call out % addValue(self % bounds % bin(i))

      ! Print upper bin boundar
      call out % addValue(self % bounds % bin(i + 1))

    end do
    call out % endArray()

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(radialMap), intent(inout) :: self

    call kill_super(self)

    call self % bounds % kill()
    self % origin = ZERO
    self % axis   = 0
    self % N      = 0

  end subroutine kill

end module radialMap_class
