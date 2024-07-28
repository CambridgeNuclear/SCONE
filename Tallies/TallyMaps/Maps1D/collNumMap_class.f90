module collNumMap_class

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
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
  interface collNumMap
    module procedure collNumMap_fromDict
  end interface

  !!
  !! Map that divides based on the number of collisions a particle underwent
  !!
  !! Private Members:
  !!   binMap  -> intMap that maps collNumber to binIdx
  !!   default -> binIdx for numbers not in binMap
  !!   Nbins   -> Number of bins in the map
  !!   collisionNumbers -> List of collision numbers required, stored in the map
  !!
  !! Interface:
  !!   tallyMap Interface
  !!   build -> builds instance without dictionary
  !!
  !! Sample Dictionary Input:
  !!   myMap {
  !!     type collNumMap;
  !!     collNumbers ( 0 1 2 3 5 10 );
  !!   }
  !!
  type, public,extends(tallyMap1D) :: collNumMap
    private
    type(intMap)                                 :: binMap
    integer(shortInt)                            :: default = 0
    integer(shortInt)                            :: Nbins   = 0
    integer(shortInt), dimension(:), allocatable :: collisionNumbers

  contains
    ! Superclass interface implementation
    procedure :: init
    procedure :: bins
    procedure :: map
    procedure :: getAxisName
    procedure :: print
    procedure :: kill

    ! Class specific procedures
    procedure :: build

  end type collNumMap

contains

  !!
  !! Build collision number map from dictionary
  !!
  !! Args:
  !!   collNumbers [in] -> Array of collision numbers to be included in the map
  !!
  !! Erorrs:
  !!   None from here.
  !!
  subroutine build(self, collNumbers)
    class(collNumMap), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: collNumbers
    integer(shortInt)                           :: i, N
    character(100), parameter :: Here = 'build (collNumMap_class.f90)'

    ! Find number of collision numbers to bin
    N = size(collNumbers)
    if (N == 0) call fatalError(Here, 'No collision number was specified')

    ! Allocate array with number of collisions
    allocate(self % collisionNumbers(N))
    self % collisionNumbers = collNumbers

    ! Allocate space in map
    call self % binMap % init(N)

    ! Load collision numbers and bins
    do i = 1, N
      call self % binMap % add(collNumbers(i), i)
    end do

    ! Save number of bins
    self % Nbins = N

  end subroutine build

  !!
  !! Initialise collision number map from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(collNumMap), intent(inout)             :: self
    class(dictionary), intent(in)                :: dict
    integer(shortInt), dimension(:), allocatable :: collNum

    ! Get cell names list
    call dict % get(collNum, 'collNumbers')

    ! Initialise Map
    call self % build(collNum)

  end subroutine init

  !!
  !! Return total number of bins in this division
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(collNumMap), intent(in) :: self
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
    class(collNumMap), intent(in)    :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    idx = self % binMap % getOrDefault( state % collisionN, self % default)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(collNumMap), intent(in) :: self
    character(nameLen)            :: name

    name = 'CollisionNumber'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(collNumMap), intent(in)    :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) // 'Bins'

    call out % startArray(name, [1, self % Nbins])

    ! Print cell indexes
    do i = 1, self % Nbins
      call out % addValue(numToChar(self % collisionNumbers(i)))
    end do

    call out % endArray()

  end subroutine print

  !!
  !! Build new collision number Map from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised collNumMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function collNumMap_fromDict(dict) result(new)
    class(dictionary), intent(in)   :: dict
    type(collNumMap)                :: new

    call new % init(dict)

  end function collNumMap_fromDict

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(collNumMap), intent(inout) :: self

    call kill_super(self)

    call self % binMap % kill()
    self % Nbins = 0

    if (allocated(self % collisionNumbers)) deallocate(self % collisionNumbers)

  end subroutine kill


end module collNumMap_class
