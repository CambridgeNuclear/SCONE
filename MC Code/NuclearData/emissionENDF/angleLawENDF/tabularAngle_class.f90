module tabularAngle_class

  use numPrecision
  use genericProcedures,  only : binarySearch, searchError, interpolate, fatalError, isSorted
  use muEndfPdf_class,    only : muEndfPdf_ptr
  use RNG_class,          only : RNG
  use angleLawENDF_class, only : angleLawENDF

  implicit none
  private

  interface tabularAngle
    module procedure new_tabularAngle
  end interface


  type,public,extends(angleLawENDF) :: tabularAngle
    !! Class to contain emission angle propabilities for secondary particles
      private
      real(defReal),dimension(:),allocatable        :: eGrid
      type(muEndfPdf_ptr),dimension(:),allocatable  :: muEndfPdfs
    contains
      procedure :: init
      procedure :: sample
      procedure :: probabilityOf
  end type tabularAngle

contains



  function sample(self,E,rand) result (mu)
    class(tabularAngle), intent(in)   :: self
    real(defReal), intent(in)         :: E
    class(RNG), intent(inout)         :: rand
    real(defReal)                     :: mu
    integer(shortInt)                 :: idx
    real(defReal)                     :: r, eps
    character(100),parameter          :: Here='sample (tabularAngle_class.f90)'

    idx = binarySearch(self % eGrid,E)
    call searchError(idx,Here)

    eps = E - self % eGrid(idx)
    r = rand % get()

    if(r < eps) then
      mu = self % muEndfPdfs(idx+1) % sample(rand)
    else
      mu = self % muEndfPdfs(idx) % sample(rand)
    end if

  end function sample


  function probabilityOf(self,mu,E) result (prob)
    class(tabularAngle), intent(in) :: self
    real(defReal), intent(in)         :: E, mu
    real(defReal)                     :: prob
    integer(shortInt)                 :: idx
    real(defReal)                     :: prob_1, prob_0, E_1, E_0
    character(100),parameter          :: Here='probabilityOf (tabularAngle_class.f90)'

    idx = binarySearch(self % eGrid,E)
    call searchError(idx,Here)

    prob_0 = self % muEndfPdfs(idx)   % probabilityOf(mu)
    prob_1 = self % muEndfPdfs(idx+1) % probabilityOf(mu)

    E_0 = self % eGrid(idx)
    E_1 = self % eGrid(idx+1)

    prob = interpolate(E_0, E_1, prob_0, prob_1, E)


  end function probabilityOf


  subroutine init(self,eGrid,muEndfPdfs)
    class(tabularAngle),intent(inout)             :: self
    real(defReal),dimension(:), intent(in)        :: eGrid
    type(muEndfPdf_ptr),dimension(:), intent(in)  :: muEndfPdfs
    character(100),parameter                      :: Here='init (tabularAngle_class.f90)'

    if(size(eGrid) /= size(muEndfPdfs)) call fatalError(Here,&
                                                        'eGrid and muEndfPdfs have diffrent size')

    if(.not.(isSorted(eGrid)))          call fatalError(Here,'eGrid is not sorted ascending')
    if(count( eGrid < 0.0 ) > 0)        call fatalError(Here,'eGrid contains -ve values')

    if(allocated(self % eGrid))  deallocate(self % eGrid)
    if(allocated(self % muEndfPdfs)) deallocate(self % muEndfPdfs)

    self % eGrid  = eGrid
    self % muEndfPdfs = muEndfPdfs

  end subroutine

  function new_tabularAngle(eGrid,muEndfPdfs) result(new)
    real(defReal),dimension(:), intent(in)        :: eGrid
    type(muEndfPdf_ptr),dimension(:), intent(in)  :: muEndfPdfs
    type(tabularAngle), pointer                   :: new

    allocate(new)
    call new % init(eGrid,muEndfPdfs)

  end function

    
end module tabularAngle_class
