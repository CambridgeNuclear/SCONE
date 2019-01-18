module tallyClerk_inter

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty

  implicit none
  private

  !!
  !! Abstract interface for a single tallyClerk.
  !! It recives reports from tallyAdmin and process them into scores on scoreMemory
  !!
  !! Its responsibilites are as follows:
  !! 1) Score some result on scoreMemory by accepting a subset of all avalible reports
  !! 2) Display implementation determined measure of convergance (usually some variance)
  !! 3) Can return information about reports it requires
  !! 4) Can return a tallyResult object for interaction with Physics Package
  !!
  !! Every tally Clerk is allocated memory location on score memory
  !! Every tally Clerk has a name
  !!
  type, public,abstract :: tallyClerk
    private
    integer(longInt)   :: memAdress = -1 ! Location of the tallyClerk on scoreMemory
    character(nameLen) :: name           ! Name of the tallyClerk

  contains
    ! Procedures used during build
    procedure(init),deferred          :: init
    procedure(validReports), deferred :: validReports

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
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: isConverged

    ! Output procedures
    procedure(display), deferred      :: display
    procedure(print),deferred         :: print
    procedure                         :: getResult

  end type tallyClerk

  !!
  !! Extandable superclass procedures
  !!
  public :: setMemAddress
  public :: setName


  abstract interface
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    function validReports(self) result(validCodes)
      import :: tallyClerk ,&
                shortInt
      class(tallyClerk),intent(in)               :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
    end function validReports

    !!
    !! Display convergance progress on the console
    !!
    subroutine display(self, mem)
      import :: tallyClerk, &
                scoreMemory
      class(tallyClerk), intent(in)    :: self
      type(scoreMemory), intent(inout) :: mem
    end subroutine display

    !!
    !! Initialise tally clerk from a dictionary
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
    subroutine print(self, outFile, mem)
      import :: tallyClerk, &
                outputFile, &
                scoreMemory
      class(tallyClerk), intent(in)    :: self
      class(outputFile), intent(inout) :: outFile
      type(scoreMemory), intent(inout) :: mem
    end subroutine print

  end interface

contains

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self,p, mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter    :: Here = 'reportInColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self, p, MT, muL, mem)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    type(scoreMemory), intent(inout)      :: mem
    character(100),parameter  :: Here = 'reportOutColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self, p, L, mem)
    class(tallyClerk), intent(inout)     :: self
    class(particle), intent(in)          :: p
    real(defReal), intent(in)            :: L
    type(scoreMemory), intent(inout)     :: mem
    character(100),parameter  :: Here = 'reportPath (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self, p, mem)
    class(tallyClerk), intent(inout)     :: self
    class(particle), intent(in)          :: p
    type(scoreMemory), intent(inout)     :: mem
    character(100),parameter  :: Here = 'reportTrans (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! Particle is associated with one of the fate codes in tallyCodes
  !!
  subroutine reportHist(self, p, mem)
    class(tallyClerk), intent(inout) :: self
    class(particle), intent(in)      :: p
    type(scoreMemory), intent(inout) :: mem
    character(100),parameter  :: Here = 'reportHist (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self, start, mem)
    class(tallyClerk), intent(inout)    :: self
    class(particleDungeon), intent(in)  :: start
    type(scoreMemory), intent(inout)    :: mem
    character(100),parameter  :: Here = 'reportCycleStart (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(tallyClerk), intent(inout)   :: self
    class(particleDungeon), intent(in) :: end
    type(scoreMemory), intent(inout)   :: mem
    character(100),parameter  :: Here = 'reportCycleEnd (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleEnd

  !!
  !! Perform convergance check in the Clerk
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
  !! Set memory adress for the clerk
  !!
  elemental subroutine setMemAddress(self, addr)
    class(tallyClerk), intent(inout) :: self
    integer(longInt), intent(in)     :: addr

    self % memAdress = addr

  end subroutine setMemAddress

  !!
  !! Return memory address of the clerk
  !!
  elemental function getMemAddress(self) result(addr)
    class(tallyClerk), intent(in) :: self
    integer(longInt)              :: addr

    addr = self % memAdress

  end function getMemAddress

  !!
  !! Return result from the clerk for interaction with Physics Package
  !! By default returns a null result
  !! Needs to be overriden in a subclass
  !!
  pure subroutine getResult(self, res, mem)
    class(tallyClerk), intent(in)                  :: self
    class(tallyResult),allocatable, intent(inout)  :: res
    type(scoreMemory), intent(inout)               :: mem

    if(allocated(res)) deallocate(res)
    allocate(tallyResultEmpty :: res)

  end subroutine getResult

  !!
  !! Set name of the clerk
  !!
  elemental subroutine setName(self, name)
    class(tallyClerk), intent(inout) :: self
    character(nameLen), intent(in)   :: name

    self % name = name

  end subroutine setName

  !!
  !! Return name of the clerk
  !!
  elemental function getName(self) result(name)
    class(tallyClerk), intent(in)  :: self
    character(nameLen)             :: name

    name = self % name

  end function getName

end module tallyClerk_inter
