module correlatedLawENDFslot_class

  use numPrecision
  use correlatedLawENDF_inter, only : correlatedLawENDF
  use RNG_class,               only : RNG

  implicit none
  private

  type, public,extends(correlatedLawENDF) :: correlatedLawENDFslot
    private
    class(correlatedLawENDF),allocatable :: slot
  contains
    ! Duplicate interface of the superclass
    procedure :: sample
    procedure :: probabilityOf

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom

  end type correlatedLawENDFslot

contains

  !!
  !! Samples mu and E_out givent incident energy E_in and random nummber generator
  !!
  subroutine sample(self,mu,E_out,E_in,rand)
    class(correlatedLawENDFslot), intent(in) :: self
    real(defReal), intent(out)               :: mu
    real(defReal), intent(out)               :: E_out
    real(defReal), intent(in)                :: E_in
    class(RNG), intent(inout)                :: rand

    call self % slot % sample(mu,E_out,E_in,rand)

  end subroutine

  !!
  !! Returns probability that neutron was emmited at mu & E_out given incident energy E_in
  !!
  function probabilityOf(self,mu,E_out,E_in) result(prob)
    class(correlatedLawENDFslot), intent(in) :: self
    real(defReal), intent(in)                :: mu
    real(defReal), intent(in)                :: E_out
    real(defReal), intent(in)                :: E_in
    real(defReal)                            :: prob

    prob = self % slot % probabilityOf(mu,E_out,E_in)

  end function probabilityOf


  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(correlatedLawENDFslot), intent(inout) :: LHS
    class(correlatedLawENDF), intent(in)        :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(correlatedLawENDFslot), intent(inout) :: LHS
    type(correlatedLawENDFslot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

    
end module correlatedLawENDFslot_class
