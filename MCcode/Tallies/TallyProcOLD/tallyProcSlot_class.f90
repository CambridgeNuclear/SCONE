module tallyProcSlot_class

  use numPrecision
  use genericProcedures,  only : fatalError
  use particle_class,     only : particle, phaseCoord
  use tallyCounter_class, only : tallyCounter
  use tallyProc_inter,    only : tallyProc

  implicit none
  private

  type, public :: tallyProcSlot
    private
    class(tallyProc),allocatable  :: proc
  contains
    ! Interface procedures
    procedure  :: filterCollision
    procedure  :: numOfBins
    procedure  :: assignBins_ptr
    procedure  :: assignBins_loc

    ! Build procedure
    generic            :: assignment(=) => assignProc
    procedure,private  :: assignProc
  end type tallyProcSlot

contains

  !!
  !! Given information about the collision and response code
  !! Return pointer to bins that need to be scored
  !!
  function filterCollision(self,pre,post,MT,muL,resCode) result (bins)
    class(tallyProcSlot), intent(in)        :: self
    class(phaseCoord), intent(in)           :: pre
    class(particle), intent(in)             :: post
    integer(shortInt), intent(in)           :: MT
    real(defReal), intent(in)               :: muL
    integer(shortInt), intent(in)           :: resCode
    type(tallyCounter),dimension(:),pointer :: bins

    bins = self % proc % filterCollision(pre,post,MT,muL,resCode)

  end function filterCollision

  !!
  !! Return number of bins tallyProc uses to store all data
  !!
  function numOfBins(self) result(N)
    class(tallyProcSlot), intent(in) :: self
    integer(shortInt)                :: N

    N = self % proc % numOfBins()

  end function numOfBins

  !!
  !! Given pointer to array of bins of size equal to tallyProc numOfBins
  !! Use the bins provided by a pointer for tallying
  !!
  subroutine assignBins_ptr(self,binMem)
    class(tallyProcSlot),intent(inout)                 :: self
    type(tallyCounter),dimension(:),pointer,intent(in) :: binMem

    call self % proc % assignBins_ptr(binMem)

  end subroutine assignBins_ptr

  !!
  !! Allocate space for tallying bins
  !!
  subroutine assignBins_loc(self)
    class(tallyProcSlot), intent(inout) :: self

    call self % proc % assignBins_loc()

  end subroutine assignBins_loc

  !!
  !! Load tallyProcessor into the slot
  !! Current content will be deallocated if present
  !! RHS will NOT be deallocated
  !!
  subroutine assignProc(LHS,RHS)
    class(tallyProcSlot), intent(inout)        :: LHS
    class(tallyProc),allocatable,intent(in)    :: RHS
    character(100),parameter                   :: Here='attachProc (tallyProcSlot_class.f90)'

    if(allocated(LHS % proc)) deallocate (LHS % proc)

    if(.not.allocated(RHS)) then
      call fatalError(Here,'Trying to load unallocated tally processor')
    end if

    !call move_alloc(RHS, LHS % proc)
    allocate(LHS % proc, source = RHS)

  end subroutine assignProc


    
end module tallyProcSlot_class
