module matXsMap_class

  use numPrecision
  use universalVariables
  use genericProcedures,       only : fatalError
  use dictionary_class,        only : dictionary
  use grid_class,              only : grid
  use particle_class,          only : particle
  use outputFile_class,        only : outputFile
  use tallyMap_inter,          only : tallyMap

  use nuclearDataRegistry_mod,    only : getMatIdx
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  implicit none
  private

  !!
  !! Constructor
  !!
  interface matXsMap
    module procedure matXsMap_fromDict
  end interface

  !!
  !! Very weird map that creates bins depending on value of total XS in a specific material
  !! Unfortunatly cannot detect upper and lower bound. Need to be provided manually.
  !! I wander if there is any practical use for it.
  !!
  type, public,extends(tallyMap) :: matXsMap
    private
    type(grid)        :: binBounds ! Bin grid
    integer(shortInt) :: N         ! Number of Bins
    integer(shortInt) :: matIdx    ! Material Index of querry material

  contains
    ! Superclass interface implementaction
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

    ! Class specific procedures
    procedure  :: init
  end type matXsMap

contains

  !!
  !! Initialise from min and max value, number of bins and extrapolation
  !!
  subroutine init(self, mini, maxi, N, type, matName)
    class(matXsMap), intent(inout) :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    character(nameLen), intent(in)  :: type
    character(nameLen), intent(in)  :: matName

    self % N = N
    call self % binBounds % init(mini, maxi, N, type)

    self % matIdx = getMatIdx(matName)

  end subroutine init


  !!
  !! Return total number of bins in this division
  !!
  pure function bins(self) result(N)
    class(matXsMap), intent(in)    :: self
    integer(shortInt)               :: N

    N = self % N

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  function map(self,p) result(idx)
    class(matXsMap), intent(in)     :: self
    class(particle), intent(in)     :: p
    integer(shortInt)               :: idx
    real(defReal)                   :: SigmaTot
    character(100),parameter :: Here = 'map (matXsMap_class.f90)'

    ! Obtain XS
    associate (xsData => p % xsData)
      select type(xsData)
        class is (transportNuclearData)
          SigmaTot = xsData % getTotalMatXS(p, self % matIdx)

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    idx = self % binBounds % search(SigmaTot)
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  function getAxisName(self) result(name)
    class(matXsMap), intent(in) :: self
    character(nameLen)              :: name

    name = 'matXs'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(matXsMap), intent(in)     :: self
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
  !! Return instance of matXsMap from dictionary
  !!
  function matXsMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(matXsMap)                         :: new
    character(nameLen)                     :: str, matName
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N
    character(100), parameter     :: Here = 'matXsMap_fromDict (matXsMap_class.f90)'

    if(.not.dict % isPresent('grid')) call fatalError(Here,"Keyword 'grid' must be present")
    if(.not.dict % isPresent('mat'))  call fatalError(Here,"Keyword 'mat' must be present")

    ! Read material name
    call dict % get(matName,'mat')



    ! Read grid definition keyword
    call dict % get(str,'grid')

    ! Choose approperiate definition
    select case(str)
      case('lin','log')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')

        ! Initialise
        call new % init(mini, maxi, N, str, matName)

      case default
        call fatalError(Here,"'grid' keyword must be: lin or log")

    end select

  end function matXsMap_fromDict

end module matXsMap_class
