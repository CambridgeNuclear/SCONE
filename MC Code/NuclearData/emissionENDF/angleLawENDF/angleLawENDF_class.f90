module angleLawENDF_class

  use numPrecision
  use genericProcedures, only : binarySearch, searchError, interpolate
  use miuEndfPdf_class,  only : miuEndfPdf_ptr
  use RNG_class,         only : RNG

  implicit none
  private

  type,public :: angleLawENDF
    !! Class to contain emission angle propabilities for secondary particles
      private
      real(defReal),dimension(:),allocatable        :: eGrid
      type(miuEndfPdf_ptr),dimension(:),allocatable :: miuEndfPdfs
    contains
      procedure :: init
      procedure :: sample
      procedure :: probabilityOf
  end type angleLawENDF

contains



  function sample(self,E,rand) result (miu)
    class(angleLawENDF), intent(in)  :: self
    real(defReal), intent(in)        :: E
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: miu
    integer(shortInt)                :: idx
    real(defReal)                    :: r, eps
    character(100),parameter         :: Here='sample (angleLawENDF_class.f90)'

    idx = binarySearch(self % eGrid,E)
    call searchError(idx,Here)

    eps = E - self % eGrid(idx)
    r = rand % get()

    if(r < eps) then
      miu = self % miuEndfPdfs(idx+1) % sample(rand)
    else
      miu = self % miuEndfPdfs(idx) % sample(rand)
    end if

  end function sample


  function probabilityOf(self,miu,E) result (prob)
    class(angleLawENDF), intent(in) :: self
    real(defReal), intent(in)       :: E, miu
    real(defReal)                   :: prob
    integer(shortInt)               :: idx
    real(defReal)                   :: prob_1, prob_0, E_1, E_0
    character(100),parameter         :: Here='probabilityOf (angleLawENDF_class.f90)'

    idx = binarySearch(self % eGrid,E)
    call searchError(idx,Here)

    prob_0 = self % miuEndfPdfs(idx)   % probabilityOf(miu)
    prob_1 = self % miuEndfPdfs(idx+1) % probabilityOf(miu)

    E_0 = self % eGrid(idx)
    E_1 = self % eGrid(idx+1)

    prob = interpolate(E_0, E_1, prob_0, prob_1, E)


  end function probabilityOf


  subroutine init(self,eGrid,miuEndfPdfs)
    class(angleLawENDF),intent(inout)             :: self
    real(defReal),dimension(:), intent(in)        :: eGrid
    type(miuEndfPdf_ptr),dimension(:), intent(in) :: miuEndfPdfs

    if(allocated(self % eGrid))  deallocate(self % eGrid)
    if(allocated(self % miuEndfPdfs)) deallocate(self % miuEndfPdfs)

    self % eGrid  = eGrid
    self % miuEndfPdfs = miuEndfPdfs
  end subroutine
    
end module angleLawENDF_class
