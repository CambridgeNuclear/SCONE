module energyMap_class

  use numPrecision
  use universalVariables
  use genericProcedures,   only : fatalError
  use dictionary_class,    only : dictionary
  use grid_class,          only : grid
  use particle_class,      only : particleState
  use outputFile_class,    only : outputFile
  use tallyMap1D_inter,    only : tallyMap1D

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
  !! Returns idx = 0 for MG particles
  !! Returns index 0 for elements outside division
  !! NOTE: Behaviour of points exactly at the boundary between two bins is undefined.
  !!       They can be in either of two bins.
  !!
  type, public,extends(tallyMap1D) :: energyMap
    private
    type(grid)        :: binBounds ! Bin grid
    integer(shortInt) :: N         ! Number of Bins

  contains
    ! Superclass interface implementaction
    procedure  :: init        ! Initialise from dictionary
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

    ! Class specific procedures
    generic            :: build => build_fromGrid, build_structured
    procedure,private  :: build_fromGrid
    procedure,private  :: build_structured
  end type energyMap

contains

  !!
  !! Build from explicit grid of bin boundaries
  !!
  subroutine build_fromGrid(self,grid)
    class(energyMap), intent(inout)         :: self
    real(defReal), dimension(:), intent(in) :: grid

    self % N = size(grid)-1
    call self % binBounds % init(grid)

  end subroutine build_fromGrid

  !!
  !! Build from min and max value, number of bins and extrapolation
  !!
  subroutine build_structured(self, mini, maxi, N, type)
    class(energyMap), intent(inout) :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    character(nameLen), intent(in)  :: type

    self % N = N
    call self % binBounds % init(mini, maxi, N, type)

  end subroutine build_structured

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(energyMap), intent(inout)        :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N
    character(100), parameter     :: Here = 'init (energyMap_class.f90)'

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
        call self % build(mini, maxi, N, type)

      case('log')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'log'

        ! Initialise
        call self % build(mini, maxi, N, type)

      case('unstruct')
        ! Read settings
        call dict % get(bins,'bins')

        ! Initialise
        call self % build(bins)

      case('predef')
        call fatalError(Here,"Predefined energy grids are not yet implemented")

      case default
        call fatalError(Here,"'grid' keyword must be: lin, log, usntruct or predef")

    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along dimension D
  !! For D=0 return all bins
  !!
  elemental function bins(self, D) result(N)
    class(energyMap), intent(in)    :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % N
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division or MG
  !!
  elemental function map(self,state) result(idx)
    class(energyMap), intent(in)     :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    ! Catch MG particle
    if( state % isMG) then
      idx = 0
      return
    end if

    ! Find position on the grid
    idx = self % binBounds % search(state % E)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  function getAxisName(self) result(name)
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

    call out % startArray(name,[self % N,2])
    do i=1,self % N
      ! Print lower bin boundary
      call out % addValue(self % binBounds % bin(i))
    end do

    do i=1,self % N
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

    call new % init(dict)

  end function energyMap_fromDict
    
end module energyMap_class
