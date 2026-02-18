module spaceMap_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use grid_class,        only : grid
  use particle_class,    only : particleState
  use outputFile_class,  only : outputFile
  use tallyMap1D_inter,  only : tallyMap1D, kill_super => kill

  implicit none
  private

  !!
  !! Constructor
  !!
  interface spaceMap
    module procedure spaceMap_fromDict
  end interface

  !!
  !! Map that divides space along one of cartesian directions (x,y or z)
  !! into number of bins
  !!
  !! Private Members:
  !!   binBounds -> grid with bin boundaries
  !!   N -> Integer number of bins in the map
  !!   dir -> direction code for the axis direction from universalVariables
  !!
  !! Interface:
  !!   tallyMap interface
  !!   build -> build instance of spaceMap without dictionary
  !!
  !! NOTE: Behaviour of points exactly at the boundary between two bins is undefined.
  !!       They can be in either of two bins.
  !!
  !! Sample Dictionary Input:
  !!   structMap {
  !!     type spaceMap;
  !!     axis x;
  !!     grid lin;
  !!     min -10.0;
  !!     max 10.0;
  !!     N 10;
  !!   }
  !!
  !!  unstructMap {
  !!    type spaceMap;
  !!    axis y;
  !!    grid unstruct;
  !!    bins (-10.0 -3.0 1.7 30.2);
  !!  }
  !!
  type, public,extends(tallyMap1D) :: spaceMap
    private
    type(grid)        :: binBounds
    integer(shortInt) :: N      = 0
    integer(shortInt) :: dir    = -17
  contains
    ! Superclass interface implementaction
    procedure  :: init
    procedure  :: bins
    procedure  :: map
    procedure  :: getAxisName
    procedure  :: print
    procedure  :: kill

    ! Class specific procedures
    generic            :: build => build_fromGrid, build_structured
    procedure,private  :: build_fromGrid
    procedure,private  :: build_structured
  end type spaceMap

contains

  !!
  !! Initialise from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(spaceMap), intent(inout)         :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N, axis
    character(100), parameter     :: Here = 'init (spaceMap_class.f90)'

    if(.not.dict % isPresent('grid')) call fatalError(Here,"Keyword 'grid' must be present")
    if(.not.dict % isPresent('axis')) call fatalError(Here,"Keyword 'axis' must be present")

    ! Find axis of division
    call dict % get(str,'axis')
    select case(str)
      case('x')
        axis = X_axis

      case('y')
        axis = Y_axis

      case('z')
        axis = Z_axis

      case default
        call fatalError(Here,'Unrecognised axis: '//trim(str)//' must be x, y or z')
    end select


    ! Read grid definition keyword
    call dict % get(str,'grid')

    ! Choose approperiate spacing definition
    select case(str)
      case('lin')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'lin'

        ! Initialise
        call self % build(mini, maxi, N, axis)

      case('unstruct')
        ! Read settings
        call dict % get(bins,'bins')

        ! Initialise
        call self % build(bins, axis)

      case default
        call fatalError(Here,"'grid' keyword must be: lin or unstruct")

    end select

  end subroutine

  !!
  !! Build from explicit grid of bin boundaries
  !!
  !! Args:
  !!   grid [in] -> Array of sorted acending defReal values
  !!   axis [in] -> Axis specification code from universalVariables
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_fromGrid(self, grid, axis)
    class(spaceMap), intent(inout)          :: self
    real(defReal), dimension(:), intent(in) :: grid
    integer(shortInt), intent(in)           :: axis

    self % N = size(grid)-1
    self % dir = axis
    call self % binBounds % init(grid)

  end subroutine build_fromGrid

  !!
  !! Build from min and max value, number of bins and direction
  !!
  !! Args:
  !!   mini [in] -> minumum value on the grid
  !!   maxi [in] -> maximum value on the grid
  !!   N [in] -> Number of bins in the grid
  !!   axis [in] -> Axis specification code from universalVariables
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_structured(self, mini, maxi, N, axis)
    class(spaceMap), intent(inout)  :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    integer(shortInt), intent(in)   :: axis
    character(nameLen)              :: type

    self % N = N
    type = 'lin'
    self % dir = axis
    call self % binBounds % init(mini, maxi, N, type)

  end subroutine build_structured

  !!
  !! Return total number of bins in this division along dimension D
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(spaceMap), intent(in)     :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % N
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification
  !!
  elemental function map(self,state) result(idx)
    class(spaceMap), intent(in)      :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx
    real(defReal)                    :: r

    r = state % r(self % dir)
    idx = self % binBounds % search(r)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(spaceMap), intent(in)     :: self
    character(nameLen)              :: name

    select case(self % dir)
      case(X_axis)
        name = 'X'
      case(Y_axis)
        name = 'Y'
      case(Z_axis)
        name ='Z'
      case default
        name ='WTF?'
    end select

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(spaceMap), intent(in)      :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bounds'

    call out % startArray(name,[2,self % N])
    do i=1,self % N
      ! Print lower bin boundary
      call out % addValue(self % binBounds % bin(i))

      ! Print upper bin boundar
      call out % addValue(self % binBounds % bin(i+1))

    end do
    call out % endArray()

  end subroutine print

  !!
  !! Return instance of spaceMap from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised spaceMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function spaceMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(spaceMap)                         :: new

    call new % init(dict)

  end function spaceMap_fromDict

  !!
  !! Kill spaceMap
  !!
  elemental subroutine kill(self)
    class(spaceMap), intent(inout) :: self

    call kill_super(self)

    call self % binBounds % kill()
    self % N = 0
    self % dir = -17

  end subroutine kill

end module spaceMap_class
