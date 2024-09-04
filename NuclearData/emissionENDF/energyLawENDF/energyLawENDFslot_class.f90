module energyLawENDFslot_class

  use numPrecision
  use energyLawENDF_inter, only : energyLawENDF
  use RNG_class,           only : RNG

  implicit none
  private

  type, public,extends(energyLawENDF) :: energyLawENDFslot
    private
    class(energyLawENDF), allocatable :: slot
  contains
    ! Duplictae interface of the superclass
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom

  end type energyLawENDFslot

contains

  !!
  !! Sample outgoing energy given random number generator and incedent energy
  !!
  function sample(self,E_in,rand) result (E_out)
    class(energyLawENDFslot), intent(in) :: self
    real(defReal), intent(in)        :: E_in
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: E_out

    E_out = self % slot % sample(E_in,rand)

  end function

  !!
  !! Give probability of outgoing energy given incedent energy
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(energyLawENDFslot), intent(in) :: self
    real(defReal), intent(in)        :: E_out,E_in
    real(defReal)                    :: prob

    prob = self % slot % probabilityOf(E_out,E_in)

  end function

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(energyLawENDFslot), intent(inout) :: self

    if(allocated(self % slot)) then
      call self % slot % kill()
      deallocate(self % slot)
    end if

  end subroutine kill

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(energyLawENDFslot), intent(inout) :: LHS
    class(energyLawENDF), intent(in)        :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(energyLawENDFslot), intent(inout)         :: LHS
    class(energyLawENDF),allocatable, intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS, LHS % slot)

  end subroutine moveAllocFrom

end module energyLawENDFslot_class
