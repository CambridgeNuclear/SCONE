module tallyClerk_inter

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty

  ! Nuclear Data Interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Extandable superclass procedures
  !!
  public :: setMemAddress
  public :: setName
  public :: kill

  !!
  !! Abstract interface for a single tallyClerk.
  !! It recives reports from tallyAdmin and process them into scores on scoreMemory
  !!
  !! Its responsibilites are as follows:
  !! 1) Score some result on scoreMemory by accepting a subset of all avalible reports
  !! 2) Display implementation determined measure of convergence (usually some variance)
  !! 3) Can return information about reports it requires
  !! 4) Can return a tallyResult object for interaction with Physics Package
  !!
  !! Every tally Clerk is allocated memory location on score memory
  !! Every tally Clerk has a name
  !!
  !! An implementation of a Clerk needs to override procedures related to reports to
  !! be able to process events reprts. Default behaviour for every report is to throw fatalError.
  !!
  !! Private Members:
  !!   memAddress -> Location of the Clerk Memory on ScoreMemory (64-bit integer)
  !!   name -> String with Clerk Name
  !!
  !! Interface:
  !!   init             -> Initialise Clerk from a dictionary
  !!   kill             -> Return to Uninitialised State
  !!   validReports     -> Returns array of integers with tallyCodes of reports that the clark requires
  !!   getSize          -> Return size required by Clerk on ScoreMemory
  !!   setMemAddress    -> Setter for "memAddress" member
  !!   getMemAddress    -> Getter for "memAddress" member
  !!   setName          -> Setter for "name" member
  !!   getName          -> Getter for "name" member
  !!   reportInColl     -> Process an incident into collision report
  !!   reportOutColl    -> Process an outgoing from collision report
  !!   reportPath       -> Process pathlength report
  !!   reportTrans      -> Process transition report
  !!   reportSpawn      -> Process particle generation report
  !!   reportHist       -> Process history report
  !!   reportCycleStart -> Process beginning of a cycle report
  !!   reportCycleEnd   -> Process end of a cycle report
  !!   closeCycle       -> Performs operations, e.g., calculate functions of scores like k-eff
  !!   isConverged      -> Return .true. if convergence criterion has been reached
  !!   display          -> Display to the console current value of a Score
  !!   print            -> Print results to the output file
  !!   getResult        -> Return tally result object for interactions
  !!
  !!
  type, public,abstract :: tallyClerk
    private
    integer(longInt)   :: memAdress = -1
    character(nameLen) :: name = ''

  contains
    ! Procedures used during build
    procedure(init),deferred          :: init
    procedure                         :: kill
    procedure(validReports), deferred :: validReports
    procedure(getSize),deferred       :: getSize

    ! Assign and get memory
    procedure                  :: setMemAddress
    procedure, non_overridable :: getMemAddress

    ! Assign an get name
    procedure                  :: setName
    procedure, non_overridable :: getName

    ! File reports and check status -> run-time procedures
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportSpawn
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: closeCycle
    procedure :: isConverged

    ! Output procedures
    procedure(display), deferred      :: display
    procedure(print),deferred         :: print
    procedure                         :: getResult

  end type tallyClerk

  abstract interface
    !!
    !! Returns array of codes that represent different reports
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Array of report tally Codes in any order (without repetitions)
    !!   Following parameters specify Report Tally Codes:
    !!     inColl_CODE
    !!     outColl_CODE
    !!     path_CODE
    !!     trans_CODE
    !!     spawn_CODE
    !!     hist_CODE
    !!     cycleStart_CODE
    !!     cycleEnd_CODE
    !!     closeCycle_CODE
    !!
    !! Errors:
    !!   None
    !!
    function validReports(self) result(validCodes)
      import :: tallyClerk ,&
                shortInt
      class(tallyClerk),intent(in)               :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
    end function validReports

    !!
    !! Return required size of score memory for the Clerk
    !!
    !! This function is crucial to calculate required size of Score Memory for the TallyAdmin
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Integer with the size of the tallyClerk
    !!
    !! Error:
    !!   Result from uninitialised Clerk is undefined
    !!
    elemental function getSize(self) result(S)
      import :: tallyClerk, &
                shortInt
      class(tallyClerk), intent(in) :: self
      integer(shortInt)             :: S
    end function getSize

    !!
    !! Display convergence progress on the console
    !!
    !! The output should aim to be kept in 60 columns
    !! If it is difficult to provide a single, sensible result
    !! it is recommended to print "<clerkClassName> does not support Display"
    !!
    !! Args:
    !!   None
    !!
    !! Errors:
    !!   None
    !!
    subroutine display(self, mem)
      import :: tallyClerk, &
                scoreMemory
      class(tallyClerk), intent(in) :: self
      type(scoreMemory), intent(in) :: mem
    end subroutine display

    !!
    !! Initialise tally clerk from a dictionary
    !!
    !! Args:
    !!   dict [in] -> Dictionary with Clerk Settings
    !!   name [in] -> Name of the Clerk
    !!
    !! Errors:
    !!   Depend on specific Clerk
    !!
    subroutine init(self, dict, name)
      import :: tallyClerk, &
                dictionary, &
                nameLen
      class(tallyClerk), intent(inout) :: self
      class(dictionary), intent(in)    :: dict
      character(nameLen), intent(in)   :: name
    end subroutine init

    !!
    !! Write contents of the clerk to output file
    !!
    !! Add results in the Clerk to the output file
    !!
    !! Args:
    !!   outFile [inout] -> Output File object
    !!   mem [in]        -> Score Memory with data
    !!
    subroutine print(self, outFile, mem)
      import :: tallyClerk, &
                outputFile, &
                scoreMemory
      class(tallyClerk), intent(in)    :: self
      class(outputFile), intent(inout) :: outFile
      type(scoreMemory), intent(in)    :: mem
    end subroutine print

  end interface

contains

  !!
  !! Process incoming collision report
  !!
  !! See tallyAdmin_class for implicit assumptions about the report.
  !!
  !! Args:
  !!   p [in]         -> Particle
  !!   xsData [inout] -> Nuclear Database with XSs data
  !!   mem [inout]    -> Score Memory to put results on
  !!   virtual [in]   -> Flag indicating virtual collision
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportInColl(self,p, xsData, mem, virtual)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual
    character(100),parameter    :: Here = 'reportInColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   p [in]        -> Particle
  !!   MT [in]       -> MT number of reaction that partilce underwent in the collision
  !!   muL [in]      -> Cosine of the collision deflection angle in LAB frame
  !!   xsData [inout]-> Nuclear Database with XSs data
  !!   mem [inout]   -> Score Memory to put results on
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportOutColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   p [in]         -> Particle
  !!   L [in]         -> Length of the path [cm]
  !!   xsData [inout] -> Nuclear Database with XSs data
  !!   mem [inout]    -> Score Memory to put results on
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportPath(self, p, L, xsData,mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    real(defReal), intent(in)             :: L
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportPath (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportPath

  !!
  !! Process transition report
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   p [in]         -> Particle
  !!   xsData [inout] -> Nuclear Database with XSs data
  !!   mem [inout]    -> Score Memory to put results on
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportTrans(self, p, xsData, mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportTrans (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportTrans

  !!
  !! Process particle creation report
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   MT [in]   -> MT number of the reaction the particle has undergone
  !!   pOld [in] -> Particle that caused the branching event
  !!   pNew [in] -> Particle state of the newly created neutron
  !!   xsData [inout] -> Nuclear Database with XSs data
  !!   mem [inout]    -> Score Memory to put results on
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportSpawn(self, MT, pOld, pNew, xsData, mem)
    class(tallyClerk), intent(inout)      :: self
    integer(shortInt), intent(in)         :: MT
    class(particle), intent(in)           :: pOld
    class(particleState), intent(in)      :: pNew
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportSpawn (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportSpawn

  !!
  !! Process history report
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   p [in]         -> Particle
  !!   xsData [inout] -> Nuclear Database with XSs data
  !!   mem [inout]    -> Score Memory to put results on
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportHist (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   end [in]    -> Particle Dungeon with particles for THIS cycle (after normalisations)
  !!   mem [inout] -> Score Memory
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportCycleStart(self, start, mem)
    class(tallyClerk), intent(inout)    :: self
    class(particleDungeon), intent(in)  :: start
    type(scoreMemory), intent(inout)    :: mem
    character(100),parameter  :: Here = 'reportCycleStart (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   end [in]    -> Particle Dungeon with particles for NEXT cycle (before any normalisation)
  !!   mem [inout] -> Score Memory
  !!
  !! Errors:
  !!   Depend on specific Clerk
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(tallyClerk), intent(inout)   :: self
    class(particleDungeon), intent(in) :: end
    type(scoreMemory), intent(inout)   :: mem
    character(100),parameter  :: Here = 'reportCycleEnd (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine reportCycleEnd

  !!
  !! Closes the cycle
  !!
  !! See tallyAdmin_class for implicit assumptionas about the report.
  !!
  !! Args:
  !!   end [in]    -> Particle Dungeon with particles for NEXT cycle (before any normalisation)
  !!   mem [inout] -> Score Memory
  !!
  !! Errors:
  !!   Depends on specific Clerk
  !!
  subroutine closeCycle(self, end, mem)
    class(tallyClerk), intent(inout)   :: self
    class(particleDungeon), intent(in) :: end
    type(scoreMemory), intent(inout)   :: mem
    character(100),parameter  :: Here = 'closeCycle (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was sent to an instance that does not support it.')

  end subroutine closeCycle

  !!
  !! Perform convergence check in the Clerk
  !!
  !! Args:
  !!   mem [in] -> Score Memory
  !!
  !! Result:
  !!   .true. if a convergence target has been reached
  !!
  !! Error:
  !!   fatalError if specific clerk does not implement this feature
  !!
  function isConverged(self, mem) result(isIt)
    class(tallyClerk), intent(in)    :: self
    type(scoreMemory), intent(inout) :: mem
    logical(defBool)                 :: isIt
    character(100),parameter  :: Here = 'isConverged (tallyClerk_inter.f90)'

    call fatalError(Here,'Convergence check is not implemented in the instance')

    ! Avoid warning
    isIt =.false.

  end function isConverged

  !!
  !! Set memory address for the clerk
  !!
  !! Args:
  !!   addr [in] -> Address in Score Memory
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine setMemAddress(self, addr)
    class(tallyClerk), intent(inout) :: self
    integer(longInt), intent(in)     :: addr

    self % memAdress = addr

  end subroutine setMemAddress

  !!
  !! Return memory address of the clerk
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Current address of the clerk in Score Memory
  !!
  !! Errors:
  !!   Result is undefined if memory address was not set yet
  !!
  elemental function getMemAddress(self) result(addr)
    class(tallyClerk), intent(in) :: self
    integer(longInt)              :: addr

    addr = self % memAdress

  end function getMemAddress

  !!
  !! Return result from the clerk for interaction with Physics Package
  !!
  !! Needs to be overriden in a subclass
  !!
  !! Args:
  !!   res [inout] -> Allocatable instance of tallyResult
  !!   mem [in]    -> Score Memory
  !!
  !! Errors:
  !!   If tallyClerk does not implement tallyResult alloctaes "res" to tallyResultEmpty class
  !!
  !!
  pure subroutine getResult(self, res, mem)
    class(tallyClerk), intent(in)                  :: self
    class(tallyResult),allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                  :: mem

    if(allocated(res)) deallocate(res)
    allocate(tallyResultEmpty :: res)

  end subroutine getResult

  !!
  !! Set name of the clerk
  !!
  !! Args:
  !!   name [in] -> Name of the Clerk
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine setName(self, name)
    class(tallyClerk), intent(inout) :: self
    character(nameLen), intent(in)   :: name

    self % name = name

  end subroutine setName

  !!
  !! Return name of the clerk
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Name of the Clerk
  !!
  !! Errors:
  !!   Undefined if name was not set yet
  !!
  elemental function getName(self) result(name)
    class(tallyClerk), intent(in)  :: self
    character(nameLen)             :: name

    name = self % name

  end function getName

  !!
  !! Return to uninitialised state
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine  kill(self)
    class(tallyClerk), intent(inout) :: self

    self % memAdress = -1
    self % name = ''

  end subroutine kill


end module tallyClerk_inter
