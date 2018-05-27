module tabularRelease_class

  use numPrecision
  use genericProcedures,    only : fatalError, isSorted
  use aceCard_class,        only : aceCard
  use endfTable_class,      only : endfTable
  use releaseLawENDF_inter, only : releaseLawENDF


  implicit none
  private


  interface tabularRelease
    module procedure new_tabularRelease_simple
    module procedure new_tabularRelease_inter
    module procedure new_tabularRelease_fromACE
  end interface

  !!
  !! Neutron release as a function of incedent energy represented as a table
  !! Multiple interpolation regions are supported
  !!
  type, public,extends(releaseLawENDF) :: tabularRelease
    private
    type(endfTable)  :: releaseTable

  contains
    generic :: init => initSimple, initInter
    procedure :: releaseAt

    procedure,private :: initSimple
    procedure,private :: initInter
  end type tabularRelease

contains

  !!
  !! Release at energy E_in
  !!
  function releaseAt(self,E_in) result(release)
    class(tabularRelease), intent(in)  :: self
    real(defReal), intent(in)          :: E_in
    real(defReal)                      :: release

    release = self % releaseTable % at(E_in)

  end function releaseAt

  !!
  !! Initialise without interpolation regions
  !! Lin-Lin single interpoaltion region is assumed
  !!
  subroutine initSimple(self, eGrid, releaseValues)
    class(tabularRelease), intent(inout)  :: self
    real(defReal),dimension(:),intent(in) :: eGrid         ! Energy Grid
    real(defReal),dimension(:),intent(in) :: releaseValues
    character(100),parameter              :: Here='initSimple (tabularRelease_class.f90)'

    ! Check if there are any -ve values in energy grid
    if ( any( eGrid < 0.0 ) ) then
      call fatalError(Here,'In the provided energy Grid some values are -ve.')
    end if

    ! Check if there are any -ve values for release
    if( any(releaseValues < 0.0 ) ) then
      call fatalError(Here,'In the provided Nu values some are -ve.')
    end if

    ! Initialise table
    call self % releaseTable % init(eGrid, releaseValues)

  end subroutine initSimple

  !!
  !! Initialise with interpolation regions
  !!
  subroutine initInter(self, eGrid, releaseValues, bounds, interENDF)
    class(tabularRelease), intent(inout)      :: self
    real(defReal),dimension(:),intent(in)     :: eGrid         ! Energy Grid
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: bounds
    integer(shortInt),dimension(:),intent(in) :: interENDF
    character(100),parameter                  :: Here='initInter (tabularRelease_class.f90)'

    ! Check if there are any -ve values in energy grid
    if ( any( eGrid < 0.0 )  ) then
      call fatalError(Here,'In the provided energy Grid some values are -ve.')
    end if

    ! Check if there are any -ve values for release
    if( any(releaseValues < 0.0 ) ) then
      call fatalError(Here,'In the provided Nu values some are -ve.')
    end if

    ! Initialise table

    call self % releaseTable % init(eGrid, releaseValues, bounds, interENDF)

  end subroutine initInter

  !!
  !! Constructor with single lin-lin interpolation region
  !!
  function new_tabularRelease_simple(eGrid, releaseValues) result(new)
    real(defReal),dimension(:),intent(in) :: eGrid
    real(defReal),dimension(:),intent(in) :: releaseValues
    type(tabularRelease),pointer          :: new

    allocate(new)
    call new % init(eGrid, releaseValues)

  end function new_tabularRelease_simple
    
  !!
  !! Constructor with multiple interpolation regions
  !!
  function new_tabularRelease_inter(eGrid, releaseValues, bounds, interENDF) result(new)
    real(defReal),dimension(:),intent(in)     :: eGrid
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: bounds
    integer(shortInt),dimension(:),intent(in) :: interENDF
    type(tabularRelease),pointer              :: new

    allocate(new)

    call new % init(eGrid, releaseValues, bounds, interENDF)

  end function new_tabularRelease_inter

  !!
  !! Constructor from ACE
  !! Head of aceCard needs to be set to the beginning of the data (KNU+1 in Table F-5 of MCNP Manual)
  !! NOTE : Defining another init for ACE would help to avoid unnecesary reallocation of memory
  !!
  function new_tabularRelease_fromACE(ACE) result(new)
    type(aceCard), intent(inout)               :: ACE
    type(tabularRelease)                       :: new
    real(defReal),dimension(:),allocatable     :: eGrid
    real(defReal),dimension(:),allocatable     :: release
    integer(shortInt),dimension(:),allocatable :: bounds
    integer(shortInt),dimension(:),allocatable :: interENDF
    integer(shortInt)                          :: NR, N
    logical(defBool)                           :: hasInterRegions

    ! Read number of interpolation regions
    NR = ACE % readInt()
    hasInterRegions = (NR /= 0)

    ! Read interpolation region data
    if( hasInterRegions) then
      bounds    = ACE % readIntArray(NR)
      interENDF = ACE % readIntArray(NR)
    end if

    ! Read rest of the data
    N       = ACE % readInt()        ! Number of energy points
    eGrid   = ACE % readRealArray(N) ! Incident energy grid
    release = ACE % readRealArray(N) ! Release values

    ! Initialise
    if( hasInterRegions) then
      call new % init(eGrid, release, bounds, interENDF)

    else
      call new % init(eGrid, release)

    end if
  end function new_tabularRelease_fromACE

end module tabularRelease_class
