module energyMap_class

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
  interface energyMap
    module procedure energyMap_fromDict
  end interface

  !!
  !! Map that divides energy into number of discrete bins
  !!
  type, public,extends(tallyMap) :: energyMap
    private
    type(grid)        :: binBounds ! Bin grid
    integer(shortInt) :: N         ! Number of Bins

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
  end type energyMap

contains

  !!
  !! Initialise from explicit grid of bin boundaries
  !!
  subroutine init_fromGrid(self,grid)
    class(energyMap), intent(inout)         :: self
    real(defReal), dimension(:), intent(in) :: grid

    self % N = size(grid)
    call self % binBounds % init(grid)

  end subroutine init_fromGrid

  !!
  !! Initialise from min and max value, number of bins and extrapolation
  !!
  subroutine init_structured(self, mini, maxi, N, type)
    class(energyMap), intent(inout) :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    character(nameLen), intent(in)  :: type

    self % N = N +1
    call self % binBounds % init(mini, maxi, N, type)

  end subroutine init_structured


  !!
  !! Return total number of bins in this division
  !!
  pure function bins(self) result(N)
    class(energyMap), intent(in)    :: self
    integer(shortInt)               :: N

    N = self % N

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,p) result(idx)
    class(energyMap), intent(in) :: self
    class(particle), intent(in)     :: p
    integer(shortInt)               :: idx

    idx = self % binBounds % search(p % E)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  pure function getAxisName(self) result(name)
    class(energyMap), intent(in) :: self
    character(nameLen)              :: name

    name = 'Energy'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(energyMap), intent(in)     :: self
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
  !! Return instance of energyMap from dictionary
  !!
  function energyMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(energyMap)                        :: new
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N
    character(100), parameter     :: Here = 'energyMap_fromDict (energyMap_class.f90)'

    if(.not.dict % isPresent('grid')) call fatalError(Here,"Keyword 'grid' must be present")

    ! Read grid definition keyword
    call dict % get(str,'grid')

    ! Choose approperiate definition
    select case(str)
      case('lin')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'lin'

        ! Initialise
        call new % init(mini, maxi, N, type)

      case('log')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'log'

        ! Initialise
        call new % init(mini, maxi, N, type)

      case('unstruct')
        ! Read settings
        call dict % get(bins,'bins')

        ! Initialise
        call new % init(bins)

      case('predef')
        call fatalError(Here,"Predefined energy grids are not yet implemented")

      case default
        call fatalError(Here,"'grid' keyword must be: lin, log, usntruct or predef")
    end select

  end function energyMap_fromDict
    
end module energyMap_class
