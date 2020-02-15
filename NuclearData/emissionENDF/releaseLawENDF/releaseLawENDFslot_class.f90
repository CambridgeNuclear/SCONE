module releaseLawENDFslot_class

  use numPrecision
  use releaseLawENDF_inter, only : releaseLawENDF
  use RNG_class,            only : RNG

  implicit none
  private

  type, public,extends(releaseLawENDF) :: releaseLawENDFslot
    private
    class(releaseLawENDF), allocatable :: slot
  contains
    ! Superclass interface
    procedure :: releaseAt
    procedure :: kill

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom

  end type releaseLawENDFslot

contains

  !!
  !! Obtain average neutron emission for incedent energy E_in
  !!
  function releaseAt(self,E_in) result(release)
    class(releaseLawENDFslot), intent(in)  :: self
    real(defReal), intent(in)              :: E_in
    real(defReal)                          :: release

    release = self % slot % releaseAt(E_in)

  end function releaseAt

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(releaseLawENDFslot), intent(inout) :: self

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
    class(releaseLawENDFslot), intent(inout) :: LHS
    class(releaseLawENDF), intent(in)        :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(releaseLawENDFslot), intent(inout) :: LHS
    type(releaseLawENDFslot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

end module releaseLawENDFslot_class
