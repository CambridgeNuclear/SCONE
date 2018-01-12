module tabularRelease_class

  use numPrecision
  use genericProcedures, only : binarySearch, endfInterpolate, interpolate,&
                                linearCeilingIdxOpen_shortInt, searchError
  use releaseLawENDF_class, only : releaseLawENDF


  implicit none
  private

  interface endfSearch
    module procedure linearCeilingIdxOpen_shortInt
  end interface

  interface tabularRelease
    module procedure newSimple_tabularRelease
    module procedure newInter_tabularRelease
  end interface


  type, public,extends(releaseLawENDF) :: tabularRelease
    private
    real(defReal), dimension(:),allocatable    :: energyPoints  !! Incoming energy grid
    real(defReal), dimension(:),allocatable    :: releaseValues !! Secondary average release points

    logical(defBool)                           :: interFlag     !! Flag to indicate presence of multiple interpolation regions

    integer(shortInt),dimension(:),allocatable :: interBound    !! Boundaries of ENDF interpolation regions
    integer(shortInt),dimension(:),allocatable :: interENDF     !! ENDF Interpolation NUmbers

  contains
    generic   :: init => initSimple, initInter
    procedure :: releaseAt

    procedure, private :: initSimple
    procedure, private :: initInter
  end type tabularRelease

contains

  subroutine initSimple(self, energyPoints, releaseValues)
    class(tabularRelease), intent(out)    :: self
    real(defReal),dimension(:),intent(in) :: energyPoints
    real(defReal),dimension(:),intent(in) :: releaseValues

    if ( allocated(self % energyPoints)) deallocate(self % energyPoints)
    if ( allocated(self % releaseValues)) deallocate(self % releaseValues)

    self % energyPoints = energyPoints
    self % releaseValues = releaseValues
    self % interFlag = .false.

  end subroutine initSimple

  subroutine initInter(self, energyPoints, releaseValues, interBound, interENDF)
    class(tabularRelease), intent(out)        :: self
    real(defReal),dimension(:),intent(in)     :: energyPoints
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: interBound
    integer(shortInt),dimension(:),intent(in) :: interENDF


    if ( allocated(self % energyPoints)) deallocate(self % energyPoints)
    if ( allocated(self % releaseValues)) deallocate(self % releaseValues)
    if ( allocated(self % interBound)) deallocate(self % interBound)
    if ( allocated(self % interENDF)) deallocate(self % interENDF)

    self % energyPoints = energyPoints
    self % releaseValues = releaseValues

    self % interFlag = .true.
    self % interBound = interBound
    self % interENDF = interENDF

  end subroutine initInter

  function releaseAt(self,energy) result(release)
    class(tabularRelease), intent(in)  :: self
    real(defReal), intent(in)          :: energy
    real(defReal)                      :: release
    integer(shortInt)                  :: bottomIdx
    integer(shortInt)                  :: interIdx
    character(100),parameter           :: Here='releaseAt (tabularRelease_class.f90)'


    bottomIdx = binarySearch(self % energyPoints, energy)
    call searchError(bottomIdx,Here)

    if (self % interFlag) then
      ! Find index of interENDF that bounds <finish this comment>
      interIdx = endfSearch(self % interBound, bottomIdx+1) ! Note the "bottomIdx+1"; topIdx is searched for.
      call searchError(interIdx,Here)

      release = endfInterpolate(self % energyPoints(bottomIdx)      , &
                                self % energyPoints(bottomIdx + 1)  , &
                                self % releaseValues(bottomIdx)     , &
                                self % releaseValues(bottomIdx + 1) , &
                                energy                              , &
                                self % interENDF(interIdx)          )
    else
      release = interpolate(self % energyPoints(bottomIdx)      , &
                            self % energyPoints(bottomIdx + 1)  , &
                            self % releaseValues(bottomIdx)     , &
                            self % releaseValues(bottomIdx + 1) , &
                            energy                              )
    endif
  end function releaseAt

  function newSimple_tabularRelease(energyPoints, releaseValues) result (new)
    real(defReal),dimension(:),intent(in) :: energyPoints
    real(defReal),dimension(:),intent(in) :: releaseValues
    type(tabularRelease),pointer          :: new

    allocate(new)
    call new % init(energyPoints, releaseValues)

  end function newSimple_tabularRelease
    
   function newInter_tabularRelease(energyPoints, releaseValues, interBound, interENDF) result (new)
    real(defReal),dimension(:),intent(in)     :: energyPoints
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: interBound
    integer(shortInt),dimension(:),intent(in) :: interENDF
    type(tabularRelease),pointer              :: new

    allocate(new)
    call new % init(energyPoints, releaseValues, interBound, interENDF)

  end function newInter_tabularRelease


end module tabularRelease_class
