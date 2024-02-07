module keffImplicitClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables
  use genericProcedures,          only : fatalError, charCmp
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
  use keffAnalogClerk_class,      only : keffResult

  ! Cache
  use ceNeutronCache_mod,         only : ceTrackingCache => trackingCache
  use mgNeutronCache_mod,         only : mgTrackingCache => trackingCache

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
  !! Private Members:
  !!   targetSTD -> Target Standard Deviation for convergance check
  !!
  !! Interface:
  !!   tallyClerk interface
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
    real(defReal)    :: targetSTD = ZERO
    ! Settings
    logical(defBool) :: virtual = .false.
  contains
    ! Duplicate interface of the tallyClerk
    ! Procedures used during build
    procedure :: init
    procedure :: kill
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
  !! See tallyClerk_inter for details
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

    ! Handle virtual collisions
    call dict % getOrDefault(self % virtual,'handleVirtual', .false.)

  end subroutine init

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(keffImplicitClerk), intent(inout) :: self

    ! Call Superclass
    call kill_super(self)

    ! Kill self
    self % targetSTD = ZERO
    self % virtual   = .false.

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(keffImplicitClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, outColl_CODE, cycleEnd_CODE, hist_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(keffImplicitClerk), intent(in) :: self
    integer(shortInt)                    :: S

    S = 5

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(keffImplicitClerk), intent(inout)  :: self
    class(particle), intent(in)              :: p
    class(nuclearDatabase),intent(inout)     :: xsData
    type(scoreMemory), intent(inout)         :: mem
    logical(defBool), intent(in)             :: virtual
    type(neutronMacroXSs)                    :: xss
    class(neutronMaterial), pointer          :: mat
    real(defReal)                            :: nuFissXS, absXS, flux
    real(defReal)                            :: s1, s2
    character(100), parameter  :: Here = 'reportInColl (keffImplicitClerk_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if (self % virtual) then

      ! Retrieve tracking cross section from cache
      ! Select over CE and MG cache, and give error if cache was not updated properly
      if (p % isMG) then
        if (mgTrackingCache(1) % G == p % G) then
          flux = p % w / mgTrackingCache(1) % xs
        else
          call fatalError(Here, 'MG tracking cache failed to update during tracking')
        end if

      else
        if (ceTrackingCache(1) % E == p % E) then
          flux = p % w / ceTrackingCache(1) % xs
        else
          call fatalError(Here, 'CE tracking cache failed to update during tracking')
        end if

      end if

    else

      if (virtual) return
      flux = p % w / xsData % getTotalMatXS(p, p % matIdx())

    end if

    ! Ensure we're not in void (could happen when scoring virtual collisions)
    if (p % matIdx() == VOID_MAT) return

    ! Get material pointer
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))
    if (.not.associated(mat)) then
      call fatalError(Here,'Unrecognised type of material was retrived from nuclearDatabase')
    end if

    ! Obtain xss
    call mat % getMacroXSs(xss, p)

    nuFissXS = xss % nuFission
    absXS    = xss % capture + xss % fission

    s1 = nuFissXS * flux
    s2 = absXS * flux

    ! Add scores to counters
    call mem % score(s1, self % getMemAddress() + IMP_PROD)
    call mem % score(s2, self % getMemAddress() + IMP_ABS)

  end subroutine reportInColl

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(keffImplicitClerk), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: MT
    real(defReal), intent(in)               :: muL
    class(nuclearDatabase),intent(inout)    :: xsData
    type(scoreMemory), intent(inout)        :: mem
    real(defReal)                           :: score

    ! Select analog score
    ! Assumes N_XNs are by implicit weight change
    select case(MT)
      case(N_2N)
        score = 1.0_defReal * p % preCollision % wgt
      case(N_3N)
        score = 2.0_defReal * p % preCollision % wgt
      case(N_4N)
        score = 3.0_defReal * p % preCollision % wgt
      case(macroAllScatter) ! Catch weight change for MG scattering
        score = max(p % w - p % preCollision % wgt, ZERO)
      case default
        score = ZERO
    end select

    ! Add to scattering production estimator
    ! Use pre collision weight
    if (score > ZERO) then
      call mem % score(score, self % getMemAddress() + SCATTER_PROD)
    end if

  end subroutine reportOutColl

  !!
  !! Process history report
  !! Gets fate code from the particle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(keffImplicitClerk), intent(inout) :: self
    class(particle), intent(in)             :: p
    class(nuclearDatabase),intent(inout)    :: xsData
    type(scoreMemory), intent(inout)        :: mem
    real(defReal)                           :: histWgt

    if( p % fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = p % w

      ! Score analog leakage
      call mem % score( histWgt, self % getMemAddress() + ANA_LEAK)

    end if

  end subroutine reportHist

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(keffImplicitClerk), intent(inout) :: self
    class(particleDungeon), intent(in)      :: end
    type(scoreMemory), intent(inout)        :: mem
    integer(longInt)                        :: addr
    real(defReal)                           :: nuFiss, absorb, leakage, scatterMul, k_est

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
  !! See tallyClerk_inter for details
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
  !! See tallyClerk_inter for details
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
  !! See tallyClerk_inter for details
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
  !! See tallyClerk_inter for details
  !!
  !! Allocates res to 'keffResult' defined in keffAnalogClerk_class
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
