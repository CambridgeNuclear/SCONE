module sphericalMap_class

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
  !! Divides space into a mesh in spherical co-ordinates
  !!
  !! TODO:
  !!   Implement polar & azimuthal subdivision
  !!
  !! Interface:
  !!   tallyMap interface
  !!
  !! NOTE:
  !!   Behaviour of points exactly at the boundary of bins is undefined.
  !!   particle can end-up in either of the two
  !!
  !! Sample Dictionary Input:
  !!  sphericalMap {
  !!    type sphericalMap;
  !!    #origin (1.0 0.0 0.0);# // Optional. Default (0 0 0)
  !!    grid lin;
  !!    #Rmin 2.0;#  // Optional. Default 0.0
  !!    Rmax  10.0;
  !!    N 10;
  !!    }
  !!
  !!
  type, public, extends (tallyMap) :: sphericalMap
    private
    real(defReal), dimension(3) :: origin = ZERO
    type(grid)                  :: rBounds
    integer(shortInt)           :: N = 0
  contains
    ! Superclass
    procedure :: init
    procedure :: bins
    procedure :: dimensions
    procedure :: getAxisName
    procedure :: map
    procedure :: print
    procedure :: kill
  end type sphericalMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(sphericalMap), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: temp, grid
    character(nameLen)                       :: type
    real(defReal)                            :: Rmin, Rmax, vol
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'init (sphericalmap_class.f90)'

    ! Check & load origin
    call dict % getOrDefault(temp, 'origin', [ZERO, ZERO, ZERO])

    if (size(temp) /= 3) then
      call fatalError(Here, 'Expected 3 values for origin. Got: ' // numToChar(size(temp)))
    end if
    self % origin = temp

    ! Load grid information
    if (.not.dict % isPresent('grid')) call fatalError(Here, 'Keyword grid must be present')
    call dict % get(type, 'grid')

    select case(type)
      case('lin')
        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % N, 'N')

        ! Check that minimum radius is OK
        if (Rmin < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve is:'//numToChar(Rmin))
        end if

        ! Build grid
        call self % rBounds % init(Rmin, Rmax, self % N, type)

      case('unstruct')

        call dict % get(grid,'bins')

        ! Initialise
        self % N = size(grid) - 1
        call self % rBounds % init(grid)

      case('equivolume')

        call dict % getOrDefault(Rmin, 'Rmin', ZERO)
        call dict % get(Rmax, 'Rmax')
        call dict % get(self % N, 'N')

        ! Check that minimum radius is OK
        if (Rmin < ZERO) then
          call fatalError(Here, 'Minumum radius must be +ve is:'//numToChar(Rmin))
        end if

        ! Calculate volume
        vol = (Rmax**3 - Rmin**3) /self % N

        allocate(grid(self % N + 1))
        ! Calculate grid boundaries
        grid(1) = Rmin
        grid(self % N + 1) = Rmax
        do i = 2,self % N
          grid(i) = (vol + grid(i-1)**3)**(ONE/3.0_defReal)
        end do

        call self % rBounds % init(grid)

      case default
        call fatalError(Here, "'grid' can take only values of: lin")

    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(sphericalMap), intent(in) :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % N
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
    class(sphericalMap), intent(in)    :: self
    integer(shortInt)                  :: D

    D = 1

  end function dimensions

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(sphericalMap), intent(in) :: self
    character(nameLen)          :: name

    name ='sphericalMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self, state) result(idx)
    class(sphericalMap), intent(in)  :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx
    real(defReal)                    :: r

    ! Calculate the distance from the origin
    r = norm2(state % r - self % origin)

    ! Search and return 0 if r is out-of-bounds
    idx = self % rBounds % search(r)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification.
  !!
  subroutine print(self,out)
    class(sphericalMap), intent(in)  :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bounds'

    call out % startArray(name, [2,self % N])
    do i = 1, self % N
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
    class(sphericalMap), intent(inout) :: self

    call kill_super(self)

    self % origin = ZERO
    call self % rBounds % kill()
    self % N = 0

  end subroutine kill

end module sphericalMap_class
