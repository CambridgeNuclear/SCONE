module equiBin32Mu_class

  use numPrecision
  use genericProcedures , only : linearFloorIdxClosed_Real, searchError
  use RNG_class , only : RNG
  use muEndfPdf_class,   only : muEndfPdf

  implicit none
  private


  interface equiBin32Mu
    module procedure new_equiBin32Mu
  end interface

 interface linSearch
    module procedure linearFloorIdxClosed_Real
  end interface


  type, public,extends(muEndfPdf) :: equiBin32Mu
    !! Class that stores PDF of mu in 32 equiprobable bins.
    private
    real(defReal),dimension(33) :: boundaries
  contains
    procedure :: sample
    procedure :: probabilityOf

    procedure,private :: init
  end type equiBin32Mu

contains

  function sample(self,rand) result (mu)
    class(equiBin32Mu), intent(in)  :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: mu
    integer(shortInt)               :: binIdx
    real(defReal)                   :: f

    binIdx = floor(32 * rand % get())
    f = rand % get()
    mu = (1.0-f)*self % boundaries(binIdx) + f* self% boundaries(binIdx+1)

  end function sample


  function probabilityOf(self,mu) result(prob)
    class(equiBin32Mu), intent(in)  :: self
    real(defReal), intent(in)       :: mu
    real(defReal)                   :: prob
    integer(shortInt)               :: binIdx
    character(100),parameter        :: Here='probabilityOf (equiBin32Mu_class.f90)'

    binIdx = linSearch(self % boundaries, mu)
    call searchError(binIdx,Here)
    prob = 1.0 / 32.0 / (self % boundaries(binIdx+1) - self % boundaries(binIdx) )

  end function probabilityOf


  subroutine init(self,boundaries)
    class(equiBin32Mu), intent(inout)        :: self
    real(defReal), dimension(33), intent(in) :: boundaries

    self % boundaries = boundaries

  end subroutine init

    
  function new_equiBin32Mu(boundaries)
    real(defReal),dimension(33), intent(in) :: boundaries
    type(equiBin32Mu),pointer               :: new_equiBin32Mu

    allocate(new_equiBin32Mu)
    call new_equiBin32Mu % init(boundaries)

  end function new_equiBin32Mu

end module equiBin32Mu_class
