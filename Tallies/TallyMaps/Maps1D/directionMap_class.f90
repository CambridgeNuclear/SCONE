module directionMap_class

  use numPrecision
  use universalVariables, only : valueOutsideArray, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, dotProduct, numToChar
  use dictionary_class,   only : dictionary
  use grid_class,         only : grid
  use particle_class,     only : particleState
  use outputFile_class,   only : outputFile
  use tallyMap1D_inter,   only : tallyMap1D, kill_super => kill

  implicit none
  private

  !!
  !! Maps the particle direction in linear bins in the range -180 -> 180 degrees (default)
  !! or boundaries defined by the user
  !!
  !! The angle is calculated using as a reference the positive direction of the
  !! axis x in the 'xy' and 'xz' cases, or y in the 'yz' case
  !!
  !! NOTE: the map is built in radians for convenience when using the function atan2, but
  !!       the user input and the results are in degrees for easy of interpretation
  !!
  !! Interface:
  !!   tallyMap1D interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   Particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  directionMap {
  !!    type directionMap;
  !!    #plane xz;#          // Optional. Default xy
  !!    N 10;
  !!    #min 60;#           // Optional. Default -180
  !!    #max 120;#          // Optional. Default 180
  !!    }
  !!
  type, public, extends (tallyMap1D) :: directionMap
    private
    type(grid)                  :: bounds
    integer(shortInt)           :: N = 0
    integer(shortInt)           :: DIM1 = 0
    integer(shortInt)           :: DIM2 = 0

  contains
    ! Superclass
    procedure :: init
    procedure :: bins
    procedure :: getAxisName
    procedure :: map
    procedure :: print
    procedure :: kill

  end type directionMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(directionMap), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(nameLen)                 :: type
    real(defReal)                      :: min, max
    character(100), parameter :: Here = 'init (directionMap_class.f90)'

    ! Check orientation of the cylinder
    if (dict % isPresent('plane')) then
      call dict % get(type, 'plane')
    else
      type = 'xy'
    end if

    select case(type)
      case('yz')
        self % DIM1 = Y_AXIS
        self % DIM2 = Z_AXIS

      case('xz')
        self % DIM1 = X_AXIS
        self % DIM2 = Z_AXIS

      case('xy')
        self % DIM1 = X_AXIS
        self % DIM2 = Y_AXIS

      case default
        call fatalError(Here, 'Keyword orientation must be xy, yz or xz. It is: '//type)

    end select

    ! Get grid details
    call dict % get(self % N, 'N')
    call dict % getOrDefault(min, 'min', -180.0_defReal)
    call dict % getOrDefault(max, 'max', 180.0_defReal)

    ! Check boundaries
    if (min < -180.0_defReal .or. min > 180.0_defReal) then
      call fatalError(Here, 'Minimum angle must be between -180 and 180 degrees. It is: '//type)
    end if

    if (max < -180.0_defReal .or. max > 180.0_defReal) then
      call fatalError(Here, 'Maximum angle must be between -180 and 180 degrees. It is: '//type)
    end if

    ! Build map in radians
    min = min*PI/180.0_defReal
    max = max*PI/180.0_defReal
    call self % bounds % init(min, max, self % N, 'lin')

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(directionMap), intent(in) :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

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
    class(directionMap), intent(in)  :: self
    character(nameLen)               :: name

    name = 'directionMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self, state) result(idx)
    class(directionMap), intent(in)  :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx
    real(defReal)                    :: x, y, theta

    ! Map the angle
    x = state % dir(self % DIM1)
    y = state % dir(self % DIM2)

    ! Returns angle in radians in the range -PI to PI
    theta = atan2(y,x)

    ! Search in the grid return 0 if index is out-of-bounds
    idx = self % bounds % search(theta)

    ! Should never happen
    if (idx == valueOutsideArray) then
      idx = 0
      return
    end if

  end function map

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification.
  !!
  subroutine print(self,out)
    class(directionMap), intent(in)  :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Print grid bins
    name = trim(self % getAxisName()) //'AngularBounds'

    call out % startArray(name, [2,self % N])

    do i = 1, self % N
      ! Print lower bin boundary in degrees
      call out % addValue(self % bounds % bin(i)/PI*180.0_defReal)

      ! Print upper bin boundary in degrees
      call out % addValue(self % bounds % bin(i + 1)/PI*180.0_defReal)
    end do

    call out % endArray()

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(directionMap), intent(inout) :: self

    call kill_super(self)

    call self % bounds % kill()
    self % N = 0
    self % DIM1 = 0
    self % DIM2 = 0

  end subroutine kill

end module directionMap_class
