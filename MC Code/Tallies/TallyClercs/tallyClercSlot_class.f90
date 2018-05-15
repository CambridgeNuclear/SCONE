module tallyClercSlot_class

  use numPrecision
!  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use tallyClerc_inter,      only : tallyClerc

  implicit none
  private


  !!
  !! Slot to store polymorphic instances of tallyClercs in an array
  !! Duplicates an interface
  !!
  type, public,extends(tallyClerc) :: tallyClercSlot
    private
    class(tallyClerc),allocatable :: slot
  contains
    ! Duplicate interface of the tallyClerc
    procedure :: reportCollision
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

  end type tallyClercSlot

contains

  !!
  !! Process collision report
  !!
  subroutine reportCollision(self,pre,post,MT,muL)
    class(tallyClercSlot), intent(inout)      :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL

    ! Pass call to instance in the slot
    call self % slot % reportCollision(pre,post,MT,muL)

  end subroutine reportCollision

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,pre,post,cellId,L)
    class(tallyClercSlot), intent(inout)     :: self
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
    class(tallyClercSlot), intent(inout)     :: self
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
    class(tallyClercSlot), intent(inout) :: self
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
    class(tallyClercSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start

    ! Pass call to instance in the slot
    call self % slot % reportCycleStart(start)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyClercSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end

    ! Pass call to instance in the slot
    call self % slot % reportCycleEnd(end)

  end subroutine reportCycleEnd

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(tallyClercSlot),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    ! Pass call to instance in the slot
    validCodes = self % slot % validReports()

  end function validReports

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self)
    class(tallyClercSlot), intent(in)  :: self

    ! Pass call to instance in the slot
    call self % slot % display()

  end subroutine display

  !!
  !! Copy alloctatable on RHS into space on LHS
  !! Does NOT deallocate RHS
  !! RHS is allocatable
  !!
  subroutine copy(LHS,RHS)
    class(tallyClercSlot), intent(inout)      :: LHS
    class(tallyClerc),allocatable,intent(in)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy
    
end module tallyClercSlot_class
