module nuclideMemoryNoMT_class

  use numPrecision
  use genericProcedures,       only : fatalError
  use xsMainSet_class,         only : xsMainSet
  use xsMainCDF_class,         only : xsMainCDF
  use xsEnergyPointNoMT_class, only : xsEnergyPointNoMT
  use aceNoMT_class,           only : aceNoMT

  implicit none
  private

  !!
  !! Extension of xsEnergyPointNoMT to store interpolated cross-section data together with
  !! key information about the nuclide. It has memory so search and interpolation is performed
  !! only when energy value requested has changed.
  !!
  type, public,extends(xsEnergyPointNoMT) :: nuclideMemoryNoMT
    private
    integer(shortInt), public        :: nucIdx       = -17     !! Nuclide index id
    real(defReal)                    :: E            = -1.0    !! Current energy of xs and pointers
    real(defReal)                    :: f            = 0.0     !! Interpolation factor for energy
    logical(defBool)                 :: isInter      = .false. !! Are XS other then total up-to-date with energy
    type(xsEnergyPointNoMT), pointer :: low          => null() !! Pointer to low energy point
    type(xsEnergyPointNoMT), pointer :: top          => null() !! Pointer to top energy point
    type(aceNoMT), pointer           :: data         => null() !! Pointer to nuclide data
  contains
    procedure :: init

    procedure :: getTotal
    procedure :: setEnergy
    procedure :: getZZid

   procedure, private :: setTo

  end type nuclideMemoryNoMT

contains


  !!
  !! Initialise nuclideMemory by connecting it to an ACE data object and providng a nuclide index
  !! to store.
  !!
  subroutine init(self,nucIdx,aceData)
    class(nuclideMemoryNoMT), intent(inout) :: self
    integer(shortInt), intent(in)           :: nucIdx
    type(aceNoMT), pointer, intent(in)      :: aceData
    character(100), parameter               :: Here = 'init (nuclideMemoryNoMT_class.f90)'

    if(.not.associated(aceData)) call fatalError(Here,'Pointer to ace data is not associated')

    self % nucIdx = nucIdx
    self % data => aceData

  end subroutine init
    
  !!
  !! Change energy to new value, calculate interpolation factor
  !!
  subroutine setTo(self,E)
    class(nuclideMemoryNoMT), intent(inout)  :: self
    real(defReal), intent(in)                :: E
    integer(shortInt)                        :: idx
    real(defReal)                            :: E_top, E_low

    ! Find index on nuclide energy grid
    idx = self % data % energyIdxFor(E)

    ! Calculate interpolation factor
    E_low = self % data % energyGrid(idx)
    E_top = self % data % energyGrid(idx+1)

    self % f = (E - E_low)/(E_top - E_low)

    ! Set pointers to the boundary energy points
    self % low => self % data % xsData(idx)
    self % top => self % data % xSData(idx+1)

    ! Set tail status flag to false
    self % isInter = .false.

  end subroutine setTo

  !!
  !! For a given energy value return total xs. Goes not interpolate the tail
  !!
  function getTotal(self,E) result (total)
    class(nuclideMemoryNoMT), intent(inout)  :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: total

    ! Check if the energy has changed
    if (self % E /= E ) then
      call self % setTo(E)
      call self % interpolateTotal(self % low, self % top, self % f)
    end if

    total = self % xs % total

  end function getTotal

  !!
  !! Interpolate all xs and cdfs in nuclide Memory to energy E
  !!
  subroutine setEnergy(self,E)
    class(nuclideMemoryNoMT), intent(inout) :: self
    real(defReal), intent(in)               :: E
    logical(defBool)                        :: tailIsNotInter
    logical(defBool)                        :: sameEnergy

    tailIsNotInter = .not.( self % isInter)
    sameEnergy     = (self % E == E)

    if (sameEnergy .and. tailIsNotInter) then
      call self % interpolateTail(self % low, self % top, self % f)

    elseif (.not. sameEnergy) then
      call self % setTo(E)
      call self % interpolate(self % low, self % top, self % f)

    end if
    ! If it is at sameEnergy with interpolated tail do nothing
    ! xs and cdf are alrady at approperiate state

  end subroutine setEnergy


  !!
  !! Return ZZid of the stored nuclide
  !!
  function getZZId(self) result(zzId)
    class(nuclideMemoryNoMT), intent(in) :: self
    character(zzIdLen)                   :: zzId

    ZZid = self % data % zzId

  end function getZZId

end module nuclideMemoryNoMT_class
