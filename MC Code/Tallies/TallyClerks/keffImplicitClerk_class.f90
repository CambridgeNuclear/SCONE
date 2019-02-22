module keffImplicitClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use genericProcedures,          only : fatalError, charCmp
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  use xsMacroSet_class,           only : xsMacroSet_ptr
  use transportNuclearData_inter, only : transportNuclearData

  use scoreMemory_class,          only : scoreMemory
  use tallyResult_class,          only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,           only : tallyClerk
  use keffAnalogClerk_class,      only : keffResult

  implicit none
  private


  !! Locations of diffrent bins wrt memory Address of the clerk
  integer(shortInt), parameter :: MEM_SIZE = 5
  integer(longInt), parameter  :: IMP_PROD     = 0 ,&  ! Implicit neutron production (from fission)
                                  SCATTER_PROD = 1 ,&  ! Analog Stattering production (N,XN)
                                  IMP_ABS      = 2 ,&  ! Implicit neutron absorbtion
                                  ANA_LEAK     = 3 ,&  ! Analog Leakage
                                  K_EFF        = 4     ! k-eff estimate

  !!
  !! A simple implicit k-eff estimator based on collison estimator of reaction rates,
  !! and an analog estimators of (N,XN) reactions and leakage
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type keffImplicitClerk;
  !!   #trigger yes/no;
  !!   #SDtarget <value>;      ! Will be read only if trigger is present and "yes"
  !!
  !! }
  !!
  type, public,extends(tallyClerk) :: keffImplicitClerk
    private
    real(defReal) :: targetSTD = ZERO
  contains
    ! Duplicate interface of the tallyClerk
    ! Procedures used during build
    procedure :: init
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportHist
    procedure :: reportCycleEnd
    procedure :: isConverged

    ! Output procedures

    procedure :: display
    procedure :: print
    procedure :: getResult
  end type keffImplicitClerk

contains

  !!
  !! Initialise from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(keffImplicitClerk), intent(inout) :: self
    class(dictionary), intent(in)           :: dict
    character(nameLen), intent(in)          :: name
    character(nameLen)                      :: chr

    ! Set name
    call self % setName(name)

    ! Configure convergance trigger
    call dict % getOrDefault(chr,'trigger','no')

    ! Read convergance target
    if( charCmp(chr,'yes')) then
      call dict % get(self % targetSTD,'SDtarget')

    end if

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(keffImplicitClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, outColl_CODE, cycleEnd_CODE, hist_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(keffImplicitClerk), intent(in) :: self
    integer(shortInt)                    :: S

    S = 5

  end function getSize

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p, mem)
    class(keffImplicitClerk), intent(inout)  :: self
    class(particle), intent(in)              :: p
    type(scoreMemory), intent(inout)         :: mem
    type(xsMacroSet_ptr)                     :: XSs
    real(defReal)                            :: totalXS, nuFissXS, absXS, flux
    real(defReal)                            :: s1, s2
    character(100), parameter  :: Here = 'reportInColl (keffActiveClerk_class.f90)'

    ! Obtain XSs
    ! Check if it dynamic type is supported
    ! If it is obtain macroscopic XSs
    ! It it isn't throw error
    select type(xsData => p % xsData)
      class is (transportNuclearData)
        call xsData % getMatMacroXS(XSs, p, p % matIdx())

      class default
        call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

    end select

    totalXS  = XSs % totalXS()
    nuFissXS = XSs % nuFissionXS()
    absXS    = XSs % captureXS() + XSs % fissionXS()

    ! Calculate flux and scores
    flux = p % w / totalXS

    s1 = nuFissXS * flux
    s2 = absXS * flux

    ! Add scores to counters
    call mem % score(s1, self % getMemAddress() + IMP_PROD)
    call mem % score(s2, self % getMemAddress() + IMP_ABS)

  end subroutine reportInColl

  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self, p, MT, muL, mem)
    class(keffImplicitClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    type(scoreMemory), intent(inout)      :: mem
    real(defReal)                         :: mult

    ! Select weight multiplier
    select case(MT)
      case(N_2N)
        mult = 1.0_defReal
      case(N_3N)
        mult = 2.0_defReal
      case(N_4N)
        mult = 3.0_defReal
      case default
        mult = ZERO
    end select

    ! Add to scattering production estimator
    ! Use pre collision weight
    if (mult /= ZERO) then
      call mem % score(p % preCollision % wgt * mult, self % getMemAddress() + SCATTER_PROD)
    end if

  end subroutine reportOutColl

  !!
  !! Process history report
  !! Gets fate code from the particle
  !!
  subroutine reportHist(self, p, mem)
    class(keffImplicitClerk), intent(inout) :: self
    class(particle), intent(in)             :: p
    type(scoreMemory), intent(inout)        :: mem
    real(defReal)                           :: histWgt

    if( p % fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = p % preHistory % wgt

      ! Score analog leakage
      call mem % score( histWgt, self % getMemAddress() + ANA_LEAK)

    end if

  end subroutine reportHist

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(keffImplicitClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    type(scoreMemory), intent(inout)     :: mem
    integer(longInt)                     :: addr
    real(defReal)                        :: nuFiss, absorb, leakage, scatterMul, k_est

    if( mem % lastCycle()) then
      addr = self % getMemAddress()
      nuFiss     = mem % getScore(addr + IMP_PROD)
      absorb     = mem % getScore(addr + IMP_ABS)
      leakage    = mem % getScore(addr + ANA_LEAK)
      scatterMul = mem % getScore(addr + SCATTER_PROD)

      k_est = nuFiss / (absorb + leakage - scatterMul )
      call mem % accumulate(k_est, addr + K_EFF)
    end if

  end subroutine reportCycleEnd

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self, mem) result(isIt)
    class(keffImplicitClerk), intent(in) :: self
    type(scoreMemory), intent(inout)     :: mem
    logical(defBool)                     :: isIt
    real(defReal)                        :: k, STD

    call mem % getResult(k, STD, self % getMemAddress() + K_EFF)

    isIt = (STD < self % targetSTD)

  end function isConverged

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(keffImplicitClerk), intent(in)  :: self
    type(scoreMemory), intent(in)         :: mem
    real(defReal)                         :: k, STD

    ! Get current k-eff estimate
    call mem % getResult(k, STD, self % getMemAddress() + K_EFF )

    ! Print to console
    print '(A,F8.5,A,F8.5)', 'k-eff (implicit): ', k, ' +/- ', STD

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(keffImplicitClerk), intent(in) :: self
    class(outputFile), intent(inout)     :: outFile
    type(scoreMemory), intent(in)        :: mem
    character(nameLen)                   :: name
    real(defReal)                        :: val, STD
    integer(longInt)                     :: addr

    call outFile % startBlock( self % getName())
    addr = self % getMemAddress()

    name = 'IMP_PROD'
    call mem % getResult(val, STD, addr + IMP_PROD)
    call outFile % printResult(val, STD, name)

    name = 'IMP_ABS'
    call mem % getResult(val, STD, addr + IMP_ABS)
    call outFile % printResult(val, STD, name)

    name = 'SCATTER_PROD'
    call mem % getResult(val, STD, addr + SCATTER_PROD)
    call outFile % printResult(val, STD, name)

    name = 'ANA_LEAK'
    call mem % getResult(val, STD, addr + ANA_LEAK)
    call outFile % printResult(val, STD, name)

    name = 'K_EFF'
    call mem % getResult(val, STD, addr + K_EFF)
    call outFile % printResult(val, STD, name)

    call outFile % endBlock()

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !!
  pure subroutine getResult(self, res, mem)
    class(keffImplicitClerk), intent(in)              :: self
    class(tallyResult), allocatable, intent(inout)    :: res
    type(scoreMemory), intent(in)                     :: mem
    real(defReal)                                     :: k, STD

    ! Get result value
    call mem % getResult(k, STD, self % getMemAddress() + K_EFF)

    allocate(res, source = keffResult([k, STD]))

  end subroutine getResult

    
end module keffImplicitClerk_class
