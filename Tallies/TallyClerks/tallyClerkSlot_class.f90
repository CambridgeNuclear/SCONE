module tallyClerkSlot_class

  use numPrecision
  use genericProcedures,      only : fatalError
  use dictionary_class,       only : dictionary
  use particle_class,         only : particle, particleState
  use particleDungeon_class,  only : particleDungeon
  use tallyClerk_inter,       only : tallyClerk, setMemAddress_super => setMemAddress, &
                                                 setName_super       => setName, &
                                                 kill_super          => kill
  use tallyClerkFactory_func, only : new_tallyClerk
  use scoreMemory_class,      only : scoreMemory
  use outputFile_class,       only : outputFile
  use tallyResult_class,      only : tallyResult

  ! Nuclear Data Interface
  use nuclearDatabase_inter,  only : nuclearDatabase

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
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! Assign and get memory
    procedure :: setMemAddress

    ! Assign an get name
    procedure :: setName

    ! File reports and check status -> run-time procedures
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportSpawn
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: isConverged

    ! Output procedures

    procedure :: display
    procedure :: print
    procedure :: getResult

    ! Subclass specific procedures
    ! Define assignment
    generic   :: assignment(=) => copy
    procedure :: copy

  end type tallyClerkSlot

contains
  !!
  !! Initialise from dictionary and name
  !! Build an instance in tallyClerkFactory and store it in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(tallyClerkSlot), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (tallyClerkSlot.f90)'

    call new_tallyClerk(self % slot, dict, name)

  end subroutine init

  !!
  !! Return to Uninitialised State
  !!
  !! Deallocates contents of the slot as well
  !!
  elemental subroutine kill(self)
    class(tallyClerkSlot), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    if(allocated(self % slot)) then
      call self % slot % kill()
      deallocate(self % slot)
    end if

  end subroutine kill


  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(tallyClerkSlot),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    ! Pass call to instance in the slot
    validCodes = self % slot % validReports()

  end function validReports

  !!
  !! Return memory size of the clerk in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(tallyClerkSlot), intent(in) :: self
    integer(shortInt)                 :: S

    S = self % slot % getSize()

  end function getSize

  !!
  !! Set memory adress for the clerk slot and stored clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental subroutine setMemAddress(self, addr)
    class(tallyClerkSLot), intent(inout) :: self
    integer(longInt), intent(in)         :: addr

    ! Call superclass procedure on itself
    call setMemAddress_super(self, addr)

    ! Call procedure in a slot
    call self % slot % setMemAddress(addr)

  end subroutine setMemAddress

  !!
  !! Set name for the clerk slot and stored clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental subroutine setName(self, name)
    class(tallyClerkSlot), intent(inout) :: self
    character(nameLen), intent(in)       :: name

    ! Call superclass procedure on itself
    call setName_super(self, name)

    ! Call procedure in a slot
    call self % slot % setName(name)

  end subroutine setName

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual

    ! Pass call to instance in the slot
    call self % slot % reportInColl(p, xsData, mem, virtual)

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem

    ! Pass call to instance in the slot
    call self % slot % reportOutColl(p, MT, muL, xsData, mem)

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportPath(self, p, L, xsData, mem)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    real(defReal), intent(in)             :: L
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem

    ! Pass call to instance in the slot
    call self % slot % reportPath(p, L, xsData, mem)

  end subroutine reportPath

  !!
  !! Process transition report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTrans(self, p, xsData, mem)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem

    ! Pass call to instance in the slot
    call self % slot % reportTrans(p, xsData, mem)

  end subroutine reportTrans

  !!
  !! Process fission report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportSpawn(self, pOld, pNew, xsData, mem)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: pOld
    class(particleState), intent(in)      :: pNew
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem

    ! Pass call to instance in the slot
    call self % slot % reportSpawn(pOld, pNew, xsData, mem)

  end subroutine reportSpawn

  !!
  !! Process history report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(tallyClerkSlot), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem

    ! Pass call to instance in the slot
    call self % slot % reportHist(p, xsData, mem)

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleStart(self, start, mem)
    class(tallyClerkSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start
    type(scoreMemory), intent(inout)     :: mem

    ! Pass call to instance in the slot
    call self % slot % reportCycleStart(start, mem)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(tallyClerkSlot), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    type(scoreMemory), intent(inout)     :: mem

    ! Pass call to instance in the slot
    call self % slot % reportCycleEnd(end, mem)

  end subroutine reportCycleEnd

  !!
  !! Perform convergance check in the Clerk
  !!
  !! See tallyClerk_inter for details
  !!
  function isConverged(self, mem) result(isIt)
    class(tallyClerkSlot), intent(in) :: self
    type(scoreMemory), intent(inout)  :: mem
    logical(defBool)                  :: isIt

    ! Pass call to instance in the slot
    isIt = self % slot % isConverged(mem)

  end function isConverged

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(tallyClerkSlot), intent(in) :: self
    type(scoreMemory), intent(in)     :: mem

    ! Pass call to instance in the slot
    call self % slot % display(mem)

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(tallyClerkSlot), intent(in) :: self
    class(outputFile), intent(inout)  :: outFile
    type(scoreMemory), intent(in)     :: mem

    call self % slot % print(outFile, mem)

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !! from the clerk in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(tallyClerkSlot), intent(in)              :: self
    class(tallyResult), allocatable, intent(inout) :: res
    type(scoreMemory), intent(in)                  :: mem

    call self % slot % getResult(res, mem)

  end subroutine getResult

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
