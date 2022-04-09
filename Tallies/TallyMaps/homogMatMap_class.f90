module homogMatMap_class

  use numPrecision
  use genericProcedures,       only : fatalError
  use dictionary_class,        only : dictionary
  use intMap_class,            only : intMap
  use particle_class,          only : particleState
  use outputFile_class,        only : outputFile
  use tallyMap1D_inter,        only : tallyMap1D, kill_super => kill

  use materialMenu_mod,        only : mm_matIdx => matIdx, mm_matName => matName

  implicit none
  private

  !!
  !! Constructor
  !!
  interface homogMatMap
    module procedure homogMatMap_fromDict
  end interface

  !!
  !! Map that divides based on the material a particle is in
  !! with the possibility of grouping some materials together
  !!
  !! NOTE: the same material can't be contained in separate bins
  !!
  !! Private Members:
  !!   binMap     -> array of intMap that maps matIdx to binIdx
  !!   matIndices -> intMap that contains unique material indeces
  !!   default    -> binIdx for materials not in binMap
  !!   Nbins      -> Number of bins in the mape
  !!   binNames   -> List of bin names for printing in output file
  !!
  !! Interface:
  !!   tallyMap Interface
  !!   build -> builds instance without dictionary
  !!
  !! Sample Dictionary Input:
  !!   myMap {
  !!     type homogMatMap;
  !!     bins (bin1 bin2);
  !!     bin1 (mat1 mat2 mat3);
  !!     bin2 (mat4 mat5);
  !!     # undefBin T; #
  !!   }
  !!
  type, public,extends(tallyMap1D) :: homogMatMap
    private
    type(intMap), dimension(:), allocatable :: binMap
    type(intMap)                            :: matIndices
    integer(shortInt)                       :: default = 0
    integer(shortInt)                       :: Nbins   = 0
    character(nameLen), dimension(:), allocatable :: binNames

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

  end type homogMatMap

contains

  !!
  !! Build homogenised material map from components directly
  !!
  !! Args:
  !!   materials [in] -> Array of material names to be included in the map
  !!   idx       [in] -> Index referring to material bin number
  !!
  !! Erorrs:
  !!   Fatal error if a material was included twice in the same bin or in two
  !!   separate bins
  !!
  subroutine build(self, materials, idx)
    class(homogMatMap), intent(inout)            :: self
    character(nameLen), dimension(:), intent(in) :: materials
    integer(shortInt), intent(in)                :: idx
    integer(shortInt)                            :: N, j, matIdx, isThere
    integer(shortInt), parameter :: PRESENT = 1
    character(100), parameter    :: Here = 'build (homogMatMap_class.f90)'

    ! Find number of materials to bin
    N = size(materials)

    ! Allocate space in map and matIndices
    call self % binMap(idx) % init(N)

    ! Load material indices and bins
    do j = 1,N

      matIdx = mm_matIdx(materials(j))
      call self % binMap(idx) % add(matIdx, PRESENT)

      ! Check if a material was included twice
      isThere = self % matIndices % getOrDefault( matIdx, self % default)
      if (isThere == 1) then
        call fatalError(Here, 'Material '//materials(j)//' was included twice')
      else
        call self % matIndices % add(matIdx, PRESENT)
      end if

    end do

  end subroutine build

  !!
  !! Initialise material map from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(homogMatMap), intent(inout)             :: self
    class(dictionary), intent(in)                 :: dict
    integer(shortInt)                             :: N, i
    character(nameLen)                            :: undefined
    character(nameLen), dimension(:), allocatable :: binNames, matNames
    character(100), parameter :: Here = 'init (homogMatMap_class.f90)'

    ! Get bin names list
    call dict % get(binNames, 'bins')

    ! Find number of materials to bin
    N = size(binNames)
    allocate(self % binMap(N), self % binNames(N))
    self % binNames = binNames

    ! Initialise map of unique material indexes
    call self % matIndices % init(1)

    ! Get setting for undefined tracking
    call dict % getOrDefault(undefined, 'undefBin', 'false')

    ! Set default and number of bins
    select case (undefined)
      case('yes','y','true','TRUE','T')
        self % Nbins   = N + 1
        self % default = N + 1

      case('no', 'n', 'false', 'FALSE', 'F')
        self % Nbins   = N
        self % default = 0
        
      case default
        call fatalError(Here, undefined//' is an unrecognised entry!')
    end select

    do i = 1,N
      ! Get material names list
      call dict % get(matNames, binNames(i))
      ! Initialise bin
      call self % build(matNames, i)

    end do

  end subroutine init

  !!
  !! Return total number of bins in this division
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(homogMatMap), intent(in)  :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % Nbins
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin. Return default value for particle out of division
  !!
  !! See tallyMap for specification
  !!
  elemental function map(self,state) result(idx)
    class(homogMatMap), intent(in)     :: self
    class(particleState), intent(in)   :: state
    integer(shortInt)                  :: idx, isThere

    do idx = 1,self % Nbins
      isThere = self % binMap(idx) % getOrDefault( state % matIdx, self % default)
      if (isThere == 1) return
    end do

    idx = self % default

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(homogMatMap), intent(in)  :: self
    character(nameLen)              :: name

    name = 'HomogMat'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(homogMatMap), intent(in)   :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bins'

    call out % startArray(name,[1,self % Nbins])

    ! Print material names
    do i = 1, size(self % binNames)
      name = self % binNames(i)
      call out % addValue(name)
    end do

    ! Print 'undefined'
    if ( self % Nbins > size(self % binMap)) then
      name = 'undefined'
      call out % addValue(name)
    end if

    call out % endArray()

  end subroutine print

  !!
  !! Build new homogenised material Map from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised homogMatMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function homogMatMap_fromDict(dict) result(new)
    class(dictionary), intent(in)               :: dict
    type(homogMatMap)                           :: new

    call new % init(dict)

  end function homogMatMap_fromDict

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(homogMatMap), intent(inout) :: self
    integer(shortInt)                 :: i

    call kill_super(self)

    do i = 1, size(self % binMap)
      call self % binMap(i) % kill()
    end do
    deallocate(self % binMap)

    self % default = 0
    self % Nbins = 0
    call self % matIndices % kill()

  end subroutine kill


end module homogMatMap_class
