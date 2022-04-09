module cellMap_class

  use numPrecision
  use genericProcedures,       only : numToChar
  use dictionary_class,        only : dictionary
  use intMap_class,            only : intMap
  use particle_class,          only : particleState
  use outputFile_class,        only : outputFile
  use tallyMap1D_inter,        only : tallyMap1D, kill_super => kill

  implicit none
  private

  !!
  !! Constructor
  !!
  interface cellMap
    module procedure cellMap_fromDict
  end interface

  !!
  !! Map that divides based on the cell a particle is in
  !!
  !! Private Members:
  !!   binMap -> intMap that maps cellIdx to binIdx
  !!   default -> binIdx for cells not in binMap
  !!   Nbins -> Number of bins in the map
  !!   cellIndices -> List of cell indices in the map
  !!
  !! Interface:
  !!   tallyMap Interface
  !!   build -> builds instance without dictionary
  !!
  !! Sample Dictionary Input:
  !!   myMap {
  !!     type cellMap;
  !!     cells (cell1 cell2);
  !!     # undefBin T; #
  !!   }
  !!
  type, public,extends(tallyMap1D) :: cellMap
    private
    type(intMap)                                  :: binMap
    integer(shortInt)                             :: default = 0
    integer(shortInt)                             :: Nbins   = 0
    integer(shortInt), dimension(:), allocatable  :: cellIndices

  contains
    ! Superclass interface implementaction
    procedure :: init
    procedure :: bins
    procedure :: map
    procedure :: getAxisName
    procedure :: print
    procedure :: kill

    ! Class specific procedures
    procedure :: build

  end type cellMap

contains

  !!
  !! Build cell map from components directly
  !!
  !! Args:
  !!   cells [in]     -> Array of cell IDs to be included in the map
  !!   trackRest [in] -> Logical. For TRUE makes extra bin for all cells not
  !!     explicitly in the map.
  !!
  !! Erorrs:
  !!   None from here. Geometry should give fatalError if a cell is given but
  !!   not defined.
  !!
  subroutine build(self, cells, trackRest)
    class(cellMap), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: cells
    logical(defBool), intent(in)                :: trackRest
    integer(shortInt)                           :: N, i

    ! Find number of materials to bin
    N = size(cells)

    ! Allocate space in map and cellIndices
    call self % binMap % init(N)
    allocate(self % cellIndices(N))

    ! Load cell indices and bins
    do i=1,N
      call self % binMap % add(cells(i), i)
      self % cellIndices(i) = cells(i)

    end do

    ! Set default and number of bins
    if(trackRest) then
      self % Nbins   = N + 1
      self % default = N + 1

    else
      self % Nbins   = N
      self % default = 0
    end if

  end subroutine build

  !!
  !! Initialise material map from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(cellMap), intent(inout)              :: self
    class(dictionary), intent(in)              :: dict
    integer(shortInt),dimension(:),allocatable :: cellIDs
    character(nameLen)                         :: undefined
    logical(defBool)                           :: trackUndefined

    ! Get material names list
    call dict % get(cellIDs, 'cells')

    ! Get setting for undefined tracking
    call dict % getOrDefault(undefined, 'undefBin', 'false')

    select case(undefined)
      case('yes','y','true','TRUE','T')
        trackUndefined = .true.

      case default
        trackUndefined = .false.
    end select

    ! Initialise Map
    call self % build(cellIDs, trackUndefined)

  end subroutine init

  !!
  !! Return total number of bins in this division
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(cellMap), intent(in)    :: self
    integer(shortInt), intent(in) :: D
    integer(shortInt)             :: N

    if (D == 1 .or. D == 0) then
      N = self % Nbins
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
    class(cellMap), intent(in)       :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    idx = self % binMap % getOrDefault( state % cellIdx, self % default)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(cellMap), intent(in)  :: self
    character(nameLen)          :: name

    name = 'Cell'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(cellMap), intent(in)       :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bins'

    call out % startArray(name,[1,self % Nbins])

    ! Print material names
    do i=1,size(self % cellIndices)
      call out % addValue(numToChar(self % cellIndices(i)))

    end do

    ! Print 'undefined'
    if ( self % Nbins > size(self % cellIndices)) then
      name = 'undefined'
      call out % addValue(name)

    end if

    call out % endArray()

  end subroutine print

  !!
  !! Build new cell Map from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised cellMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function cellMap_fromDict(dict) result(new)
    class(dictionary), intent(in)   :: dict
    type(cellMap)                   :: new

    call new % init(dict)

  end function cellMap_fromDict

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cellMap), intent(inout) :: self

    call kill_super(self)

    call self % binMap % kill()
    self % default = 0
    self % Nbins = 0

  end subroutine kill


end module cellMap_class
