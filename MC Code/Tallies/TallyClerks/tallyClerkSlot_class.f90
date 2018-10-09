module tallyClerkSlot_class

  use numPrecision
  use genericProcedures,     only : fatalError
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use tallyClerk_inter,      only : tallyClerk
  use outputFile_class,      only : outputFile

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
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    procedure :: isConverged

    procedure :: validReports
    procedure :: display
    ! Initialisation procedure -> return error Slot is a container
    procedure :: init
    procedure :: print

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
  subroutine reportOutColl(self,p,MT,muL)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL

    ! Pass call to instance in the slot
    call self % slot % reportOutColl(p,MT,muL)

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,p,L)
    class(tallyClerkSlot), intent(inout) :: self
    class(particle), intent(in)          :: p
    real(defReal), intent(in)            :: L

    ! Pass call to instance in the slot
    call self % slot % reportPath(p,L)

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,p)
    class(tallyClerkSlot), intent(inout) :: self
    class(particle), intent(in)          :: p

    ! Pass call to instance in the slot
    call self % slot % reportTrans(p)

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,p)
    class(tallyClerkSlot), intent(inout) :: self
    class(particle), intent(in)          :: p

    ! Pass call to instance in the slot
    call self % slot % reportHist(p)

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
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(tallyClerkSlot), intent(in) :: self
    logical(defBool)                  :: isIt

    ! Pass call to instance in the slot
    isIt = self % slot % isConverged()

  end function isConverged

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
  !! Initialise from dictionary.
  !! Return error.
  !!
  subroutine init(self,dict,name)
    class(tallyClerkSlot), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (tallyClerkSlot.f90)'

    call fatalError(Here,'tallyClerkSlot cannot be initialised from dictionary')

  end subroutine init

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  subroutine print(self,outFile)
    class(tallyClerkSlot), intent(in) :: self
    class(outputFile), intent(inout)  :: outFile

    call self % slot % print(outFile)

  end subroutine print

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(tallyClerkSlot), intent(inout)      :: LHS
    class(tallyClerk),intent(in)              :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    allocate(LHS % slot, source = RHS)

  end subroutine copy
    
end module tallyClerkSlot_class
