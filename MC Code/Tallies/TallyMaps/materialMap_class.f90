module materialMap_class

  use numPrecision
  use genericProcedures,       only : fatalError
  use dictionary_class,        only : dictionary
  use intMap_class,            only : intMap
  use particle_class,          only : particleState
  use outputFile_class,        only : outputFile
  use tallyMap1D_inter,        only : tallyMap1D


  use nuclearDataRegistry_mod, only : getMatIdx, getMatName

  implicit none
  private

  !!
  !! Constructor
  !!
  interface materialMap
    module procedure materialMap_fromDict
  end interface

  !!
  !! Map that divides based on the material a particle is in
  !!
  type, public,extends(tallyMap1D) :: materialMap
    private
    type(intMap)                                  :: binMap
    integer(shortInt)                             :: default = 0
    integer(shortInt)                             :: Nbins
    integer(shortInt), dimension(:), allocatable  :: matIndices

  contains
    ! Superclass interface implementaction
    procedure :: init        !! Initialise from dictionary
    procedure :: bins        !! Return number of bins
    procedure :: map         !! Map particle to a bin
    procedure :: getAxisName !! Return character describing variable of devision
    procedure :: print       !! Print values associated with bins to outputfile

    ! Class specific procedures
    procedure :: build       !! Build from components (without dictionary)

  end type materialMap

contains

  !!
  !! Build material map from components directly
  !!
  subroutine build(self, materials, trackRest)
    class(materialMap), intent(inout)            :: self
    character(nameLen), dimension(:), intent(in) :: materials
    logical(defBool),intent(in)                  :: trackRest
    integer(shortInt)                            :: N, i, matIdx

    ! Find number of materials to bin
    N = size(materials)

    ! Allocate space in map and matIndices
    call self % binMap % init(N)
    allocate(self % matIndices(N))

    ! Load material indices and bins
    do i=1,N
      matIdx = getMatIdx(materials(i))
      call self % binMap % add(matIdx, i)
      self % matIndices(i) = matIdx

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
  subroutine init(self, dict)
    class(materialMap), intent(inout)           :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen),dimension(:),allocatable :: matNames
    character(nameLen)                          :: undefined
    logical(defBool)                            :: trackUndefined

    ! Get material names list
    call dict % get(matNames, 'materials')

    ! Get setting for undefined tracking
    call dict % getOrDefault(undefined, 'undefBin', 'false')

    select case(undefined)
      case('yes','y','true','TRUE','T')
        trackUndefined = .true.

      case default
        trackUndefined = .false.
    end select

    ! Initialise Map
    call self % build(matNames, trackUndefined)

  end subroutine init

  !!
  !! Return total number of bins in this division
  !!
  elemental function bins(self, D) result(N)
    class(materialMap), intent(in)  :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % Nbins
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,state) result(idx)
    class(materialMap), intent(in)     :: self
    class(particleState), intent(in)   :: state
    integer(shortInt)                  :: idx

    idx = self % binMap % getOrDefault( state % matIdx, self % default)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  function getAxisName(self) result(name)
    class(materialMap), intent(in)  :: self
    character(nameLen)              :: name

    name = 'Material'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(materialMap), intent(in)   :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bins'

    call out % startArray(name,[1,self % Nbins])

    ! Print material names
    do i=1,size(self % matIndices)
      name = getMatName(self % matIndices(i))
      call out % addValue(name)

    end do

    ! Print 'undefined'
    if ( self % Nbins > size(self % matIndices)) then
      name = 'undefined'
      call out % addValue(name)

    end if

    call out % endArray()

  end subroutine print

  !!
  !! Build new material Map from dictionary
  !!
  function materialMap_fromDict(dict) result(new)
    class(dictionary), intent(in)               :: dict
    type(materialMap)                           :: new

    call new % init(dict)

  end function materialMap_fromDict

end module materialMap_class
