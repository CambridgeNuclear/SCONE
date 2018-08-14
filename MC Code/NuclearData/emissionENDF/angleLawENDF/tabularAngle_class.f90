module tabularAngle_class

  use numPrecision
  use genericProcedures,   only : binarySearch, searchError, interpolate, fatalError, isSorted
  use aceCard_class,       only : aceCard
  use RNG_class,           only : RNG
  use angleLawENDF_inter,  only : angleLawENDF

  ! Diffrent mu pdfs
  use muEndfPdf_inter,     only : muEndfPdf
  use muEndfPdfSlot_class, only : muEndfPdfSlot
  use isotropicMu_class,   only : isotropicMu
  use equiBin32Mu_class,   only : equiBin32Mu
  use tabularMu_class,     only : tabularMu

  implicit none
  private

  interface tabularAngle
    module procedure new_tabularAngle
    module procedure new_tabularAngle_fromACE
  end interface

  !!
  !! Contains energy dependant mu data
  !!
  type,public,extends(angleLawENDF) :: tabularAngle
      private
      real(defReal),dimension(:),allocatable        :: eGrid
      type(muEndfPdfSlot),dimension(:),allocatable  :: muEndfPdfs
    contains
      procedure :: init
      procedure :: sample
      procedure :: probabilityOf
  end type tabularAngle

contains

  !!
  !! Given collison energy and random number generator sample mu
  !!
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

    eps = (E - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))
    r = rand % get()

    if(r < eps) then
      mu = self % muEndfPdfs(idx+1) % sample(rand)
    else
      mu = self % muEndfPdfs(idx) % sample(rand)

    end if

  end function sample

  !!
  !! Return probability density of mu at collision energy E
  !!
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

  !!
  !! Initialise from energy grid and array of corresponding mu PDFs at single energy
  !! NOTE : Content in muEndfPdfs slots will be deallocated (moved allocation)
  !!
  subroutine init(self,eGrid,muEndfPdfs)
    class(tabularAngle),intent(inout)                :: self
    real(defReal),dimension(:), intent(in)           :: eGrid
    type(muEndfPdfSlot),dimension(:), intent(inout)  :: muEndfPdfs
    integer(shortInt)                                 :: i
    character(100),parameter                    :: Here='init (tabularAngle_class.f90)'

    ! Perform checks
    if(size(eGrid) /= size(muEndfPdfs)) call fatalError(Here,'eGrid and muEndfPdfs have diffrent size')

    if(.not.(isSorted(eGrid)))    call fatalError(Here,'eGrid is not sorted ascending')
    if(any( eGrid < 0.0 ))        call fatalError(Here,'eGrid contains -ve values')

    if(allocated(self % eGrid))      deallocate(self % eGrid)
    if(allocated(self % muEndfPdfs)) deallocate(self % muEndfPdfs)

    ! Copy energy grid and move allocation of muEndfPdfSlots
    self % eGrid  = eGrid

    allocate(self % muEndfPdfs(size(eGrid)))

    do i=1,size(muEndfPdfs)
      call self % muEndfPdfs(i) % moveAllocFrom( muEndfPdfs(i) )

    end do

  end subroutine init

  !!
  !! Constructor of tabularAngle
  !! NOTE : Content in muEndfPdfs slots will be deallocated (moved allocation)
  !!
  function new_tabularAngle(eGrid,muEndfPdfs) result(new)
    real(defReal),dimension(:), intent(in)          :: eGrid
    type(muEndfPdfSlot),dimension(:), intent(inout) :: muEndfPdfs
    type(tabularAngle)                              :: new

    call new % init(eGrid,muEndfPdfs)

  end function new_tabularAngle

  !!
  !! Constructoe of tabularAngle from ACE
  !! ACE head should be set to beegining of tabular mu data
  !!
  function new_tabularAngle_fromACE(ACE) result(new)
    type(aceCard), intent(inout)                   :: ACE
    type(tabularAngle)                             :: new
    real(defReal), dimension(:), allocatable       :: eGrid
    integer(shortInt)                              :: N, i
    integer(shortInt),dimension(:),allocatable     :: muLoc
    type(muEndfPdfSlot),dimension(:), allocatable  :: muPdfs
    
    ! Read initial information
    N     = ACE % readInt()        ! Read size of the energy grid
    eGrid = ACE % readRealArray(N) ! Read energy grid
    muLoc = ACE % readIntArray(N)  ! Read mu pdf locators

    ! Allocate space for angleLawENDFslots
    allocate(muPdfs(N))

    ! Build array of muPdfs
    do i=1,size(muLoc)
      select case(muLoc(i))
        case(0) ! Isotropic mu pdf
          muPdfs(i) = isotropicMu()

        case(1:) ! +ve -> 32 equiprobable bin distribution
          call ACE % setToAnglePdf(muLoc(i))
          muPdfs(i) = equiBin32Mu(ACE)

        case(:-1) ! -ve -> tabular pdf
          call ACE % setToAnglePdf( abs(muLoc(i)))
          muPdfs(i) = tabularMu(ACE) ! Explicity use CDF in ACE data

        case default ! Clearly this should never happen. But codes surprise you...
          call fatalError('new_tabularAngle_fromACE (tabularAngle_class.f90)','Impossible state. WTF?')

      end select
    end do

    ! Initialise new tabularAngle
    call new % init(eGrid, muPdfs)

  end function new_tabularAngle_fromACE

end module tabularAngle_class
