module equiBin32Miu_class

  use numPrecision
  use genericProcedures , only : linearFloorIdxClosed_Real, searchError
  use RNG_class , only : RNG
  use miuEndfPdf_class,   only : miuEndfPdf

  implicit none
  private


  interface equiBin32Miu
    module procedure new_equiBin32Miu
  end interface

 interface linSearch
    module procedure linearFloorIdxClosed_Real
  end interface


  type, public,extends(miuEndfPdf) :: equiBin32Miu
    !! Class that stores PDF of miu in 32 equiprobable bins.
    private
    real(defReal),dimension(33) :: boundaries
  contains
    procedure :: sample
    procedure :: probabilityOf

    procedure,private :: init
  end type equiBin32Miu

contains

  function sample(self,rand) result (miu)
    class(equiBin32Miu), intent(in) :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: miu
    integer(shortInt)               :: binIdx
    real(defReal)                   :: f

    binIdx = floor(32 * rand % get())
    f = rand % get()
    miu = (1.0-f)*self % boundaries(binIdx) + f* self% boundaries(binIdx+1)

  end function sample


  function probabilityOf(self,miu) result(prob)
    class(equiBin32Miu), intent(in) :: self
    real(defReal), intent(in)       :: miu
    real(defReal)                   :: prob
    integer(shortInt)               :: binIdx
    character(100),parameter        :: Here='probabilityOf (equiBin32Miu_class.f90)'

    binIdx = linSearch(self % boundaries, miu)
    call searchError(binIdx,Here)
    prob = 1.0 / 32.0 / (self % boundaries(binIdx+1) - self % boundaries(binIdx) )

  end function probabilityOf


  subroutine init(self,boundaries)
    class(equiBin32Miu), intent(inout)       :: self
    real(defReal), dimension(33), intent(in) :: boundaries

    self % boundaries = boundaries

  end subroutine init

    
  function new_equiBin32Miu(boundaries)
    real(defReal),dimension(33), intent(in) :: boundaries
    type(equiBin32Miu),pointer              :: new_equiBin32Miu

    allocate(new_equiBin32Miu)
    call new_equiBin32Miu % init(boundaries)

  end function new_equiBin32Miu

end module equiBin32Miu_class
