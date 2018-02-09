module contTabularEnergy_class

  use numPrecision
  use genericProcedures,   only : binarySearch, fatalError, interpolate, searchError, isSorted
  use tabularEnergy_class, only : tabularEnergy
  use RNG_class,           only : RNG
  use energyLawENDF_class, only : energyLawENDF


  implicit none
  private

  interface contTabularEnergy
    module procedure new_contTabularEnergy
  end interface

  type, public,extends(energyLawENDF):: contTabularEnergy
    private
    real(defReal),dimension(:),allocatable       :: eGrid
    type(tabularEnergy),dimension(:),allocatable :: ePdfs
  contains
    procedure :: sample
    procedure :: probabilityOf

    procedure :: init

  end type contTabularEnergy

contains

  function sample(self,E_in,rand) result (E_out)
    class(contTabularEnergy), intent(in) :: self
    real(defReal), intent(in)            :: E_in
    class(RNG), intent(inout)            :: rand
    real(defReal)                        :: E_out
    integer(shortInt)                    :: idx
    real(defReal)                        :: r, eps
    character(100),parameter             :: Here='sample (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    eps = E_in - self % eGrid(idx)
    r = rand % get()

    if(r < eps) then
      E_out = self % ePdfs(idx+1) % sample(rand)
    else
      E_out = self % ePdfs(idx) % sample(rand)
    end if

  end function sample


  function probabilityOf(self,E_out,E_in) result (prob)
    class(contTabularEnergy), intent(in) :: self
    real(defReal), intent(in)            :: E_out, E_in
    real(defReal)                        :: prob
    integer(shortInt)                    :: idx
    real(defReal)                        :: prob_1, prob_0, E_1, E_0
    character(100),parameter             :: Here='probabilityOf (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    prob_0 = self % ePdfs(idx)   % probabilityOf(E_out)
    prob_1 = self % ePdfs(idx+1) % probabilityOf(E_out)

    E_0 = self % eGrid(idx)
    E_1 = self % eGrid(idx+1)

    prob = interpolate(E_0, E_1, prob_0, prob_1, E_in)

  end function probabilityOf


  subroutine init(self,eGrid,ePdfs)
    class(contTabularEnergy), intent(inout)     :: self
    real(defReal),dimension(:),intent(in)       :: eGrid
    type(tabularEnergy),dimension(:),intent(in) :: ePdfs
    character(100),parameter                    :: Here='init (contTabularEnergy_class.f90)'

    ! Check if the provided eGrid and ePdfs match in size and if eGrid is sorted and all its
    ! elements are +ve.
    if(size(eGrid) /= size(ePdfs))  call fatalError(Here,'eGrid and ePdfs have diffrent size')
    if(.not.(isSorted(eGrid)))      call fatalError(Here,'eGrid is not sorted ascending')
    if ( count( eGrid < 0.0 ) > 0 ) call fatalError(Here,'eGrid contains -ve values')


    if(allocated(self % eGrid)) deallocate(self % eGrid)
    if(allocated(self % ePdfs)) deallocate(self % ePdfs)

    self % eGrid  = eGrid
    self % ePdfs = ePdfs


  end subroutine init


  function new_contTabularEnergy(eGrid,ePdfs) result (new)
    real(defReal),dimension(:),intent(in)        :: eGrid
    type(tabularEnergy),dimension(:),intent(in)  :: ePdfs
    type(contTabularEnergy),pointer              :: new

    allocate(new)
    call new % init(eGrid,ePdfs)

  end function new_contTabularEnergy

end module contTabularEnergy_class
