module spaceMap_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use grid_class,        only : grid
  use particle_class,    only : particle
  use outputFile_class,  only : outputFile
  use tallyMap_inter,    only : tallyMap

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
  type, public,extends(tallyMap) :: spaceMap
    private
    type(grid)        :: binBounds ! Bin grid
    integer(shortInt) :: N         ! Number of Bins
    integer(shortInt) :: dir =-17  ! Direction code
  contains
    ! Superclass interface implementaction
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

    ! Class specific procedures
    generic            :: init => init_fromGrid, init_structured
    procedure,private  :: init_fromGrid
    procedure,private  :: init_structured
  end type spaceMap

contains

  !!
  !! Initialise from explicit grid of bin boundaries
  !!
  subroutine init_fromGrid(self,grid,axis)
    class(spaceMap), intent(inout)          :: self
    real(defReal), dimension(:), intent(in) :: grid
    integer(shortInt), intent(in)           :: axis

    self % N = size(grid)
    self % dir = axis
    call self % binBounds % init(grid)

  end subroutine init_fromGrid

  !!
  !! Initialise from min and max value, number of bins and extrapolation
  !!
  subroutine init_structured(self, mini, maxi, N, axis)
    class(spaceMap), intent(inout)  :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    integer(shortInt), intent(in)   :: axis
    character(nameLen)              :: type

    self % N = N +1
    type = 'lin'
    self % dir = axis
    call self % binBounds % init(mini, maxi, N, type)

  end subroutine init_structured


  !!
  !! Return total number of bins in this division
  !!
  pure function bins(self) result(N)
    class(spaceMap), intent(in)     :: self
    integer(shortInt)               :: N

    N = self % N

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,p) result(idx)
    class(spaceMap), intent(in)     :: self
    class(particle), intent(in)     :: p
    integer(shortInt)               :: idx

    associate (r => p % rGlobal())
      idx = self % binBounds % search( r(self % dir) )
      if (idx == valueOutsideArray) idx = 0
    end associate

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  pure function getAxisName(self) result(name)
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
  function spaceMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(spaceMap)                         :: new
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N, axis
    character(100), parameter     :: Here = 'spaceMap_fromDict (spaceMap_class.f90)'

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
        call new % init(mini, maxi, N, axis)

      case('unstruct')
        ! Read settings
        call dict % get(bins,'bins')

        ! Initialise
        call new % init(bins, axis)

      case default
        call fatalError(Here,"'grid' keyword must be: lin or usntruct")

    end select

  end function spaceMap_fromDict
    
end module spaceMap_class
