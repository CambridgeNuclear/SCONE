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
  !! Maps the particle direction in linear bins in the range 0-360 degrees
  !!
  !! Interface:
  !!   tallyMap1D interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  directionMap {
  !!    type directionMap;
  !!    #plane xz;#            // Optional. Default xy
  !!    N 10;
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
    class(directionMap), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    character(nameLen)                       :: type
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
        call fatalError(Here, 'Keyword orientation must be x, y or z. It is: '//type)

    end select

    ! Build grid
    call dict % get(self % N, 'N')
    call self % bounds % init(-PI, PI, self % N, 'lin')

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

    x = state % dir(self % DIM1)
    y = state % dir(self % DIM2)

    theta = atan2(y,x)
    ! Search along the azimuthal dimension and return 0 if index is out-of-bounds
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
    class(directionMap), intent(inout) :: self

    call kill_super(self)

    call self % bounds % kill()
    self % N = 0
    self % DIM1 = 0
    self % DIM2 = 0

  end subroutine kill

end module directionMap_class
