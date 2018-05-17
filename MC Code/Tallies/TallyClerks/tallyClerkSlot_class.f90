module tallyClerkSlot_class

  use numPrecision
!  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use tallyClerk_inter,      only : tallyClerk

  implicit none
  private


  !!
  !! Slot to store polymorphic instances of tallyClerks in an array
  !! Duplicates an interface
  !!
  type, public,extends(tallyClerk) :: tallyClerkSlot
    private
    class(tallyClerk),allocatable :: slot
  contains
    ! Duplicate interface of the tallyClerk
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    procedure :: validReports
    procedure :: display
    !procedure :: print *** Interface for this procedure will be defined shortly

    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy

  end type tallyClerkSlot

contains

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self,p)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p

    ! Pass call to instance in the slot
    call self % slot % reportInColl(p)

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self,pre,post,MT,muL)
    class(tallyClerkSlot), intent(inout)      :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL

    ! Pass call to instance in the slot
    call self % slot % reportOutColl(pre,post,MT,muL)

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,pre,post,cellId,L)
    class(tallyClerkSlot), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt), intent(in)        :: cellId
    real(defReal), intent(in)            :: L

    ! Pass call to instance in the slot
    call self % slot % reportPath(pre,post,cellId,L)

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,pre,post)
    class(tallyClerkSlot), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post

    ! Pass call to instance in the slot
    call self % slot % reportTrans(pre,post)

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,pre,post,fate)
    class(tallyClerkSlot), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate

    ! Pass call to instance in the slot
    call self % slot % reportHist(pre,post,fate)

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyClerkSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start

    ! Pass call to instance in the slot
    call self % slot % reportCycleStart(start)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyClerkSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end

    ! Pass call to instance in the slot
    call self % slot % reportCycleEnd(end)

  end subroutine reportCycleEnd

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(tallyClerkSlot),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    ! Pass call to instance in the slot
    validCodes = self % slot % validReports()

  end function validReports

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self)
    class(tallyClerkSlot), intent(in)  :: self

    ! Pass call to instance in the slot
    call self % slot % display()

  end subroutine display

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loasing slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(tallyClerkSlot), intent(inout)      :: LHS
    class(tallyClerk),intent(in)              :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy
    
end module tallyClerkSlot_class
