module removalTimeClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables,         only : VOID_MAT, TRACKING_XS, MAX_COL
  use genericProcedures,          only : fatalError, charCmp
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Nuclear Data Interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast
  use neutronXSPackages_class,    only : neutronMacroXSs

  ! Tally Interfaces
  use scoreMemory_class,          only : scoreMemory
  use tallyResult_class,          only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  implicit none
  private


  !! Locations of different bins wrt memory Address of the clerk
  integer(shortInt), parameter :: MEM_SIZE = 4
  integer(longInt), parameter  :: IMP_POP  = 0 ,&  ! Implicit population estimator
                                  IMP_ABS  = 1 ,&  ! Implicit absorption estimator
                                  ANA_LEAK = 2 ,&  ! Analog Leakage
                                  TAU      = 3     ! Removal time estimate
  !!
  !! A simple implicit removal time (tau) estimator based on collison estimator of reaction rates,
  !! and on analog estimators of leakage
  !!
  !! tau = population / (Absorption Rate + Leakage Rate)
  !!
  !! Private Members:
  !!   targetSTD -> Target Standard Deviation for convergence check
  !!
  !! Interface:
  !!   tallyClerk interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type removalTimeClerk;
  !!   #trigger yes/no;
  !!   #SDtarget <value>;      ! Will be read only if trigger is present and "yes"
  !!
  !! }
  !!
  type, public,extends(tallyClerk) :: removalTimeClerk
    private
    real(defReal)    :: targetSTD = ZERO
    ! Settings
    logical(defBool) :: handleVirtual = .true.
  contains
    ! Duplicate interface of the tallyClerk
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: closeCycle
    procedure :: isConverged

    ! Output procedures
    procedure :: display
    procedure :: print

  end type removalTimeClerk

contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(removalTimeClerk), intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen), intent(in)         :: name
    character(nameLen)                     :: chr

    ! Set name
    call self % setName(name)

    ! Configure convergence trigger
    call dict % getOrDefault(chr,'trigger','no')

    ! Read convergence target
    if( charCmp(chr,'yes')) then
      call dict % get(self % targetSTD,'SDtarget')

    end if

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .true.)

  end subroutine init

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(removalTimeClerk), intent(inout) :: self

    ! Call Superclass
    call kill_super(self)

    ! Kill self
    self % targetSTD = ZERO
    self % handleVirtual = .true.

  end subroutine kill

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(removalTimeClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, hist_CODE, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(removalTimeClerk), intent(in) :: self
    integer(shortInt)                   :: S

    S = 4

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(removalTimeClerk), intent(inout)  :: self
    class(particle), intent(in)             :: p
    class(nuclearDatabase),intent(inout)    :: xsData
    type(scoreMemory), intent(inout)        :: mem
    logical(defBool), intent(in)            :: virtual
    type(neutronMacroXSs)                   :: xss
    class(neutronMaterial), pointer         :: mat
    real(defReal)                           :: invSpeed, absXS, flux
    real(defReal)                           :: s1, s2
    character(100), parameter  :: Here = 'reportInColl (removalTimeClerk_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! Calculate flux with the right cross section according to virtual collision handling
    if (self % handleVirtual) then
      flux = p % w / xsData % getTrackingXS(p, p % matIdx(), TRACKING_XS)
    else
      flux = p % w / xsData % getTotalMatXS(p, p % matIdx())
    end if

    ! Get material pointer
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))
    if (.not.associated(mat)) then
      call fatalError(Here,'Unrecognised type of material was retrived from nuclearDatabase')
    end if

    ! Obtain xss
    call mat % getMacroXSs(xss, p)

    invSpeed = ONE / p % getSpeed()
    absXS    = xss % capture + xss % fission

    s1 = invSpeed * flux
    s2 = absXS * flux

    ! Add scores to counters
    call mem % score(s1, self % getMemAddress() + IMP_POP)
    call mem % score(s2, self % getMemAddress() + IMP_ABS)

  end subroutine reportInColl

  !!
  !! Process history report
  !! Gets fate code from the particle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(removalTimeClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    class(nuclearDatabase),intent(inout)   :: xsData
    type(scoreMemory), intent(inout)       :: mem
    real(defReal)                          :: histWgt

    if (p % fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = p % w

      ! Score analog leakage
      call mem % score(histWgt, self % getMemAddress() + ANA_LEAK)

    end if

  end subroutine reportHist

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine closeCycle(self, end, mem)
    class(removalTimeClerk), intent(inout) :: self
    class(particleDungeon), intent(in)     :: end
    type(scoreMemory), intent(inout)       :: mem
    integer(longInt)                       :: addr
    real(defReal)                          :: pop, absorb, leakage, tau_est

    if (mem % lastCycle()) then

      addr = self % getMemAddress()
      pop        = mem % getScore(addr + IMP_POP)
      absorb     = mem % getScore(addr + IMP_ABS)
      leakage    = mem % getScore(addr + ANA_LEAK)

      tau_est = pop / (absorb + leakage)
      call mem % accumulate(tau_est, addr + TAU)

    end if

  end subroutine closeCycle

  !!
  !! Perform convergence check in the Clerk
  !!
  !! See tallyClerk_inter for details
  !!
  function isConverged(self, mem) result(isIt)
    class(removalTimeClerk), intent(in) :: self
    type(scoreMemory), intent(inout)    :: mem
    logical(defBool)                    :: isIt
    real(defReal)                       :: t, STD

    call mem % getResult(t, STD, self % getMemAddress() + TAU)

    isIt = (STD < self % targetSTD)

  end function isConverged

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(removalTimeClerk), intent(in) :: self
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: t, STD
    character(MAX_COL)                  :: buffer

    ! Get current tau estimate
    call mem % getResult(t, STD, self % getMemAddress() + TAU)

    ! Print to console
    write(buffer, '(A,ES12.5,A,ES12.5)') 'Removal time: ', t, ' s +/- ', STD
    call statusMsg(buffer)

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(removalTimeClerk), intent(in) :: self
    class(outputFile), intent(inout)    :: outFile
    type(scoreMemory), intent(in)       :: mem
    character(nameLen)                  :: name
    real(defReal)                       :: val, STD
    integer(longInt)                    :: addr

    call outFile % startBlock( self % getName())
    addr = self % getMemAddress()

    name = 'IMP_POP'
    call mem % getResult(val, STD, addr + IMP_POP)
    call outFile % printResult(val, STD, name)

    name = 'IMP_ABS'
    call mem % getResult(val, STD, addr + IMP_ABS)
    call outFile % printResult(val, STD, name)

    name = 'ANA_LEAK'
    call mem % getResult(val, STD, addr + ANA_LEAK)
    call outFile % printResult(val, STD, name)

    name = 'TAU'
    call mem % getResult(val, STD, addr + TAU)
    call outFile % printResult(val, STD, name)

    call outFile % endBlock()

  end subroutine print

end module removalTimeClerk_class
