module muEndfPdfSlot_class

  use numPrecision
  use RNG_class,       only : RNG
  use muEndfPdf_inter, only : muEndfPdf

  implicit none
  private

  type, public,extends(muEndfPdf) :: muEndfPdfSlot
    private
    class(muEndfPdf),allocatable :: slot
  contains
    ! Duplicate interface of the superclass
    procedure :: sample
    procedure :: probabilityOf

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom

  end type muEndfPdfSlot

contains

  !!
  !! Sample mu given random number generator
  !!
  function sample(self,rand) result(mu)
    class(muEndfPdfSlot),intent(in)  :: self
    class(RNG),intent(inout)         :: rand
    real(defReal)                    :: mu

    mu = self % slot % sample(rand)

  end function sample

  !!
  !! Give probability density of mu
  !! Does not check if mu is in <-1;1>
  !!
  function probabilityOf(self,mu) result(prob)
    class(muEndfPdfSlot), intent(in)  :: self
    real(defReal), intent(in)         :: mu
    real(defReal)                     :: prob

    prob = self % slot % probabilityOf(mu)

  end function probabilityOf

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(muEndfPdfSlot), intent(inout) :: LHS
    class(muEndfPdf), intent(in)        :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(muEndfPdfSlot), intent(inout) :: LHS
    type(muEndfPdfSlot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

end module muEndfPdfSlot_class
