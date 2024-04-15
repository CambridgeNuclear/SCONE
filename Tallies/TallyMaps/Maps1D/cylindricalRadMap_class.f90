module cylindricalRadMap_class

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
  !! Divides space into a radial mesh in cylindrical co-ordinates
  !!
  !! Interface:
  !!   tallyMap1D interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  cylindricalRadMap {
  !!    type cylindricalRadMap;
  !!    #orientation x;#    // Optional. Default z
  !!    #origin (1.0 0.0);# // Optional. Default (0 0). Order is xy, xz or yz.
  !!    rGrid lin;
  !!    #Rmin 2.0;#         // Optional. Default 0.0
  !!    Rmax  10.0;
  !!    rN 10;
  !!    }
  !!
  type, public, extends (tallyMap1D) :: cylindricalRadMap
    private
    real(defReal), dimension(2) :: origin = ZERO
    type(grid)                  :: rBounds
    integer(shortInt)           :: rN   = 0
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

  end type cylindricalRadMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(cylindricalRadMap), intent(inout)  :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: temp, grid
    character(nameLen)                       :: type
    real(defReal)                            :: Rmin, Rmax, vol
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'init (cylindricalRadMap_class.f90)'

    ! Check & load origin
    call dict % getOrDefault(temp, 'origin', [ZERO, ZERO])

    if (size(temp) /= 2) then
      call fatalError(Here, 'Expected 2 values for origin. Got: ' // numToChar(size(temp)))
    end if
    self % origin = temp

    ! Check orientation of the cylinder
    if (dict % isPresent('orientation')) then
      call dict % get(type, 'orientation')
    else
      type = 'z'
    end if

    select case(type)
      case('x')
        self % DIM1 = Y_AXIS
        self % DIM2 = Z_AXIS

      case('y')
        self % DIM1 = X_AXIS
        self % DIM2 = Z_AXIS

      case('z')
        self % DIM1 = X_AXIS
        self % DIM2 = Y_AXIS

      case default
        call fatalError(Here, 'Keyword orientation must be x, y or z. It is: '//type)

      end select

    ! Load radial grid information
    if (.not. dict % isPresent('rGrid')) call fatalError(Here, 'Keyword rGrid must be present')
    call dict % get(type, 'rGrid')

    ! Check type of radial grid bins
    select case(type)
      case('lin')
        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % rN, 'rN')

        ! Check that minimum radius is OK
        if (Rmin < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve. It is: '//numToChar(Rmin))
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
          call fatalError(Here, 'Minumum radius must be +ve. It is: '//numToChar(Rmin))
        end if

        ! Calculate volume
        vol = (Rmax**2 - Rmin**2) /self % rN

        allocate(grid(self % rN + 1))
        ! Calculate grid boundaries
        grid(1) = Rmin
        grid(self % rN + 1) = Rmax
        do i = 2,self % rN
          grid(i) = sqrt(vol + grid(i-1)**2)
        end do

        call self % rBounds % init(grid)

      case default
        call fatalError(Here, "'rGrid' can take only values of: lin, unstruct, equivolume")

    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(cylindricalRadMap), intent(in) :: self
    integer(shortInt), intent(in)        :: D
    integer(shortInt)                    :: N

    if (D == 1 .or. D == 0) then
      N = self % rN
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
    class(cylindricalRadMap), intent(in) :: self
    character(nameLen)                   :: name

    name = 'cylindricalRadMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self, state) result(idx)
    class(cylindricalRadMap), intent(in) :: self
    class(particleState), intent(in)     :: state
    integer(shortInt)                    :: idx
    real(defReal)                        :: r

    ! Calculate the distance from the origin
    r = norm2(state % r([self % DIM1, self % DIM2]) - self % origin)

    ! Search and return 0 if r is out-of-bounds
    idx = self % rBounds % search(r)

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
    class(cylindricalRadMap), intent(in) :: self
    class(outputFile), intent(inout)     :: out
    character(nameLen)                   :: name
    integer(shortInt)                    :: i

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

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cylindricalRadMap), intent(inout) :: self

    call kill_super(self)

    self % origin = ZERO
    call self % rBounds % kill()
    self % rN   = 0
    self % DIM1 = 0
    self % DIM2 = 0

  end subroutine kill

end module cylindricalRadMap_class
