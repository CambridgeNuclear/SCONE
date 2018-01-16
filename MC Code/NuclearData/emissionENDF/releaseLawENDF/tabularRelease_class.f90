module tabularRelease_class

  use numPrecision
  use endfTable_class, only : endfTable
  use releaseLawENDF_class, only : releaseLawENDF


  implicit none
  private


  interface tabularRelease
    module procedure newSimple_tabularRelease
    module procedure newInter_tabularRelease
  end interface


  type, public,extends(releaseLawENDF) :: tabularRelease
    private
    class(endfTable),pointer  :: releaseTable => null()

  contains
    generic :: init => initSimple, initInter
    procedure :: releaseAt

    procedure,private :: initSimple
    procedure,private :: initInter
  end type tabularRelease

contains

  function releaseAt(self,energy) result(release)
    class(tabularRelease), intent(in)  :: self
    real(defReal), intent(in)          :: energy
    real(defReal)                      :: release

    release = self % releaseTable % at(energy)

  end function releaseAt


  subroutine initSimple(self, eGrid, releaseValues)
    class(tabularRelease), intent(inout)  :: self
    real(defReal),dimension(:),intent(in) :: eGrid         ! Energy Grid
    real(defReal),dimension(:),intent(in) :: releaseValues

    if(associated(self % releaseTable)) deallocate(self % releaseTable)

    self % releaseTable => endfTable(eGrid, releaseValues)

  end subroutine initSimple

   subroutine initInter(self, eGrid, releaseValues, bounds, interENDF)
    class(tabularRelease), intent(inout)      :: self
    real(defReal),dimension(:),intent(in)     :: eGrid         ! Energy Grid
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: bounds
    integer(shortInt),dimension(:),intent(in) :: interENDF

    if(associated(self % releaseTable)) deallocate(self % releaseTable)

    self % releaseTable => endfTable(eGrid, releaseValues, bounds, interENDF)

  end subroutine initInter


  function newSimple_tabularRelease(eGrid, releaseValues) result (new)
    real(defReal),dimension(:),intent(in) :: eGrid
    real(defReal),dimension(:),intent(in) :: releaseValues
    type(tabularRelease),pointer          :: new

    allocate(new)
    call new % init(eGrid, releaseValues)

  end function newSimple_tabularRelease
    

  function newInter_tabularRelease(eGrid, releaseValues, bounds, interENDF) result (new)
    real(defReal),dimension(:),intent(in)     :: eGrid
    real(defReal),dimension(:),intent(in)     :: releaseValues
    integer(shortInt),dimension(:),intent(in) :: bounds
    integer(shortInt),dimension(:),intent(in) :: interENDF
    type(tabularRelease),pointer              :: new

    allocate(new)

    call new % init(eGrid, releaseValues, bounds, interENDF)

  end function newInter_tabularRelease


end module tabularRelease_class
