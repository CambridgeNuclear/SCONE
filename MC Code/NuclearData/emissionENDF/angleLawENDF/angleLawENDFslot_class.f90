module angleLawENDFslot_class

  use numPrecision
  use angleLawENDF_inter, only : angleLawENDF
  use RNG_class,          only : RNG

  implicit none
  private

  type, public,extends(angleLawENDF) :: angleLawENDFslot
    private
    class(angleLawENDF),allocatable :: slot
  contains
    ! Duplicate interfacte of the superclass
    procedure :: sample
    procedure :: probabilityOf

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom

  end type angleLawENDFslot

contains

  !!
  !! Given collison energy and random number generator sample mu
  !!
  function sample(self,E,rand) result (mu)
    class(angleLawENDFslot), intent(in) :: self
    real(defReal), intent(in)           :: E
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: mu

    mu = self % slot % sample(E,rand)

  end function sample

  !!
  !! Return probability density of mu at collision energy E
  !!
  function probabilityOf(self,mu,E) result (prob)
    class(angleLawENDFslot), intent(in) :: self
    real(defReal), intent(in)           :: E, mu
    real(defReal)                       :: prob

    prob = self % slot % probabilityOf(mu,E)

  end function probabilityOf

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(angleLawENDFslot), intent(inout) :: LHS
    class(angleLawENDF), intent(in)        :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(angleLawENDFslot), intent(inout) :: LHS
    type(angleLawENDFslot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

end module angleLawENDFslot_class
