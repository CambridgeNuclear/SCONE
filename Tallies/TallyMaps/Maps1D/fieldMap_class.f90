module fieldMap_class

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use dictionary_class,        only : dictionary
  use particle_class,          only : particleState
  use coord_class,             only : coordList
  use outputFile_class,        only : outputFile
  use tallyMap1D_inter,        only : tallyMap1D, kill_super => kill

  ! Field
  use pieceConstantField_inter,       only : pieceConstantField
  use pieceConstantFieldFactory_func, only : new_pieceConstantField

  implicit none
  private
  
  !!
  !! Constructor
  !!
  interface fieldMap
    module procedure fieldMap_fromDict
  end interface

  !!
  !! Map that wraps fields to locate a particle. Performs the field mapping,
  !! returning the index that a particle occupies.
  !!
  !! Private Members:
  !!   default -> binIdx for cells not in binMap
  !!   Nbins   -> Number of bins in the map
  !!
  !! Interface:
  !!   tallyMap Interface
  !!
  !! Sample Dictionary Input:
  !!   myMap {
  !!     type fieldMap;
  !!     <piece constant field dictionary>
  !!
  !!   }
  !!
  type, public,extends(tallyMap1D) :: fieldMap
    private
    integer(shortInt)                      :: default = 0
    integer(shortInt)                      :: Nbins   = 0
    class(pieceConstantField), allocatable :: field

  contains
    ! Superclass interface implementaction
    procedure :: init
    procedure :: bins
    procedure :: map
    procedure :: getAxisName
    procedure :: print
    procedure :: kill

  end type fieldMap

contains

  !!
  !! Initialise cell map from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(fieldMap), intent(inout) :: self
    class(dictionary), intent(in)  :: dict
    class(dictionary),pointer      :: tempDict

    tempDict => dict % getDictPtr('field')
    call new_pieceConstantField(self % field, tempDict)
    self % Nbins = self % field % getSize()

  end subroutine init

  !!
  !! Return total number of bins in the map/elements in the field
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(fieldMap), intent(in)   :: self
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
    class(fieldMap), intent(in)      :: self
    class(particleState), intent(in) :: state
    type(coordList)                  :: coords
    integer(shortInt)                :: idx

    call coords % assignPosition(state % r)
    idx = self % field % map(coords)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(fieldMap), intent(in) :: self
    character(nameLen)          :: name

    name = 'Field'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(fieldMap), intent(in)      :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) // 'Bins'

    call out % startArray(name, [1, self % Nbins])

    ! Print field indexes
    do i = 1, self % Nbins
      call out % addValue(numToChar(i))
    end do

    call out % endArray()

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fieldMap), intent(inout) :: self

    call kill_super(self)
    call self % field % kill()
    self % Nbins = 0

  end subroutine kill
  
  !!
  !! Build new field Map from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised fieldMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function fieldMap_fromDict(dict) result(new)
    class(dictionary), intent(in) :: dict
    type(fieldMap)                :: new

    call new % init(dict)

  end function fieldMap_fromDict

end module fieldMap_class
