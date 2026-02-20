module energyMap_class

  use numPrecision
  use universalVariables
  use preDefEnergyGrids
  use genericProcedures,   only : fatalError
  use dictionary_class,    only : dictionary
  use grid_class,          only : grid
  use particle_class,      only : particleState
  use outputFile_class,    only : outputFile
  use tallyMap1D_inter,    only : tallyMap1D, kill_super => kill

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
  !!
  !! Private Members:
  !!   binBounds -> grid with bin boundaries
  !!   N -> Integer number of bins in the map
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
  !!     type energyMap;
  !!     grid log;
  !!     min 1.0E-9;
  !!     max 20.0;
  !!     N 10;
  !!   }
  !!
  !!   unstructMap {
  !!     type energyMap;
  !!     grid ustruct;
  !!     bins (1.0E-9 1.0E-8 0.6E-6 0.3 20.0);
  !!   }
  !!
  !!   predefMap {
  !!     type energyMap;
  !!     grid predef;
  !!     name wims69;
  !!   }
  !!
  type, public,extends(tallyMap1D) :: energyMap
    private
    type(grid)        :: binBounds
    integer(shortInt) :: N = 0

  contains
    ! Superclass interface implementaction
    procedure  :: init
    procedure  :: bins
    procedure  :: map
    procedure  :: getAxisName
    procedure  :: print
    procedure  :: kill

    ! Class specific procedures
    generic            :: build => build_fromGrid, build_structured, build_predef
    procedure, private :: build_fromGrid
    procedure, private :: build_structured
    procedure, private :: build_predef
    procedure          :: printReverse

  end type energyMap

contains

  !!
  !! Build from explicit grid of bin boundaries
  !!
  !! Args:
  !!   grid [in] -> Array of sorted ascending defReal values
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_fromGrid(self, grid)
    class(energyMap), intent(inout)         :: self
    real(defReal), dimension(:), intent(in) :: grid

    self % N = size(grid)-1
    call self % binBounds % init(grid)

  end subroutine build_fromGrid

  !!
  !! Build from min and max value, number of bins and direction
  !!
  !! Args:
  !!   mini [in] -> minumum value on the grid
  !!   maxi [in] -> maximum value on the grid
  !!   N [in] -> Number of bins in the grid
  !!   type [in] -> nameLen type of interpolation 'lin' or 'log'
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
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
  !! Retrieves predefined bin boundaries, reverts the order of the grid to be
  !! from thermal to fast, and calls the function to build from bins
  !!
  !! Args:
  !!   name [in] -> name of the bin boundaries array, as present in preDefEnergyGrids
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_predef(self, name)
    class(energyMap), intent(inout)        :: self
    character(nameLen), intent(in)         :: name
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: i, j, N
    real(defReal)                          :: temp
    character(100), parameter :: Here = 'build_predef (energyMap_class.f90)'

    select case(name)
      case('wims69')
        bins = wims69

      case('wims172')
        bins = wims172

      case('casmo40')
        bins = casmo40

      case('casmo23')
        bins = casmo23

      case('casmo12')
        bins = casmo12

      case('casmo7')
        bins = casmo7

      case('ecco33')
        bins = ecco33

      case('vitaminj')
        bins = vitaminj

      case default
          call fatalError(Here,'Grid '//name//' is undefined!')
    end select

    ! Revert order from low to high (thermal to fast)
    N = size(bins)
    do i = 1, N/2       ! Head moves forward
      j = N - i + 1     ! j = Tail moves backward
      temp  = bins(i)      ! swap elements
      bins(i)  = bins(j)
      bins(j)  = temp
    end do

    ! Initialise
    call self % build(bins)

  end subroutine build_predef

  !!
  !! Initialise from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(energyMap), intent(inout)        :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N
    character(nameLen)                     :: name
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
        ! Read settings
        call dict % get(name,'name')

        ! Initialise
        call self % build(name)

      case default
        call fatalError(Here,"'grid' keyword must be: lin, log, usntruct or predef")

    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along dimension D
  !! For D=0 return all bins
  !!
  !! See tallyMap for specification
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
  !! See tallyMap for specification
  !!
  !! NOTE:
  !!   Returns idx = 0 for MG particles
  !!
  elemental function map(self,state) result(idx)
    class(energyMap), intent(in)     :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    ! Catch MG particle
    if (state % isMG) then
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
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(energyMap), intent(in) :: self
    character(nameLen)              :: name

    name = 'Energy'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self, out)
    class(energyMap), intent(in)     :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bounds'

    call out % startArray(name,[self % N,2])
    do i = 1,self % N
      ! Print lower bin boundary
      call out % addValue(self % binBounds % bin(i))
    end do

    do i = 1,self % N
      ! Print upper bin boundar
      call out % addValue(self % binBounds % bin(i+1))
    end do

    call out % endArray()

  end subroutine print

  !!
  !! Add information about division axis to the output file
  !!
  !! Prints the map in the opposite order: from fast to thermal
  !!
  subroutine printReverse(self, out)
    class(energyMap), intent(in)     :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Protect from trying to print in an uninitialised state
    if (self % N == 0) return

    ! Name the array
    name = trim(self % getAxisName()) //'Bounds'

    call out % startArray(name,[self % N,2])
    do i = 1, self % N
      ! Print lower bin boundary
      call out % addValue(self % binBounds % bin(self % N - i + 1))
    end do

    do i = 1, self % N
      ! Print upper bin boundar
      call out % addValue(self % binBounds % bin(self % N - i + 2))
    end do

    call out % endArray()

  end subroutine printReverse

  !!
  !! Return instance of energyMap from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised energyMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function energyMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(energyMap)                        :: new

    call new % init(dict)

  end function energyMap_fromDict

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(energyMap), intent(inout) :: self

    call kill_super(self)

    ! Kill local
    call self % binBounds % kill()
    self % N = 0

  end subroutine kill

end module energyMap_class
