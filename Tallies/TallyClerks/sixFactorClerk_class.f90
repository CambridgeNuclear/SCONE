module sixFactorClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables,         only : VOID_MAT, TRACKING_XS, MAX_COL
  use genericProcedures,          only : fatalError, charCmp, numToChar
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


  !! Default thermal energy cutoff
  real(defReal), parameter :: DEFAULT_THERMAL = 0.6E-6

  !! Locations of different bins wrt memory Address of the clerk
  integer(shortInt), parameter :: MEM_SIZE = 14
  integer(longInt), parameter, public  :: FAST_PROD = 0 ,&
                                          THERMAL_PROD = 1, &
                                          FAST_LEAK = 2, &
                                          FAST_ABS = 3 ,&
                                          THERMAL_LEAK = 4, &
                                          THERMAL_ABS = 5 ,&
                                          THERMAL_FUEL_ABS = 6, &
                                          FAST_FISSION_FACTOR = 7 ,&
                                          FAST_NON_LEAKAGE = 8 ,&
                                          RESONANCE_ESCAPE_PROB = 9 ,&
                                          THERMAL_NON_LEAKAGE = 10 ,&
                                          THERMAL_UTILISATION = 11 ,&
                                          ETA_THERMAL = 12, &
                                          K_EFFECTIVE = 13

  !!
  !! A simple implicit six factor formula estimator based on implicit estimation
  !! of capture and fission rates, and analog estimation of (N,XN) and leakage.
  !! (N,XN) is included in eta or the fast fission factor.
  !!
  !! Private Members:
  !!   thermalEnergy -> energy boundary between fast and thermal
  !!
  !! Interface:
  !!   tallyClerk interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type sixFactorClerk;
  !!   #thermalEnergy 0.6E-6; # ! Sets the E boundary between fast and thermal
  !! }
  !!
  type, public,extends(tallyClerk) :: sixFactorClerk
    private
    ! Settings
    logical(defBool) :: handleVirtual = .true.
    real(defReal)    :: thermalEnergy = DEFAULT_THERMAL
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
    procedure :: closeCycle

    ! Output procedures
    procedure :: display
    procedure :: print

  end type sixFactorClerk

contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(sixFactorClerk), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100), parameter  :: Here = 'init (sixFactorClerk_class.f90)'

    ! Set name
    call self % setName(name)

    ! Set thermal/fast boundary
    call dict % getOrDefault(self % thermalEnergy, 'thermalEnergy', DEFAULT_THERMAL)
    if (self % thermalEnergy <= 0) then
      call fatalError(Here,'Must have a positive thermal energy boundary.')
    end if

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .true.)

  end subroutine init

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(sixFactorClerk), intent(inout) :: self

    ! Call Superclass
    call kill_super(self)

    ! Kill self
    self % thermalEnergy = DEFAULT_THERMAL
    self % handleVirtual = .true.

  end subroutine kill

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(sixFactorClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, outColl_CODE, hist_CODE, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(sixFactorClerk), intent(in) :: self
    integer(shortInt)                    :: S

    S = MEM_SIZE

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(sixFactorClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase),intent(inout)  :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual
    type(neutronMacroXSs)                 :: xss
    class(neutronMaterial), pointer       :: mat
    real(defReal)                         :: nuFissXS, absXS, flux
    real(defReal)                         :: s1, s2
    integer(shortInt)                     :: absIdx, fissProdIdx
    character(100), parameter  :: Here = 'reportInColl (sixFactorClerk_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! Ensure we're not in void (could happen when scoring virtual collisions)
    if (p % matIdx() == VOID_MAT) return

    ! Cannot currently be used for MG calculations
    if (p % isMG) return

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
    nuFissXS = xss % nuFission
    absXS = (xss % capture + xss % fission)
    s1 = nuFissXS * flux
    s2 = absXS * flux

    ! Determine whether fast or thermal and tally the appropriate quantities
    if (p % E <= self % thermalEnergy) then
      fissProdIdx = THERMAL_PROD
      absIdx = THERMAL_ABS
      
      if (mat % isFissile()) then
        call mem % score(s2, self % getMemAddress() + THERMAL_FUEL_ABS)
      end if

    else
      fissProdIdx = FAST_PROD
      absIdx = FAST_ABS
    end if

    ! Add scores to counters
    call mem % score(s1, self % getMemAddress() + fissProdIdx)
    call mem % score(s2, self % getMemAddress() + absIdx)

  end subroutine reportInColl

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(sixFactorClerk), intent(inout) :: self
    class(particle), intent(in)          :: p
    integer(shortInt), intent(in)        :: MT
    real(defReal), intent(in)            :: muL
    class(nuclearDatabase),intent(inout) :: xsData
    type(scoreMemory), intent(inout)     :: mem
    real(defReal)                        :: score

    ! Select analog score
    ! Assumes N_XNs are by implicit weight change
    select case(MT)
      case(N_2N, N_2Nd, N_2Na, N_2N2a, N_2Np, N_2Nl(1):N_2Ncont)
        score = 1.0_defReal * p % preCollision % wgt
      case(N_3N, N_3Na, N_3Np)
        score = 2.0_defReal * p % preCollision % wgt
      case(N_4N)
        score = 3.0_defReal * p % preCollision % wgt
      case(macroAllScatter, macroIEScatter)
        score = max(p % w - p % preCollision % wgt, ZERO)
      case default
        score = ZERO
    end select

    ! Add to fast fission production
    ! Use pre collision weight
    if (score > ZERO) then
      call mem % score(score, self % getMemAddress() + FAST_PROD)
    end if

  end subroutine reportOutColl

  !!
  !! Process history report
  !! Gets fate code from the particle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(sixFactorClerk), intent(inout) :: self
    class(particle), intent(in)          :: p
    class(nuclearDatabase),intent(inout) :: xsData
    type(scoreMemory), intent(inout)     :: mem
    real(defReal)                        :: histWgt
    integer(shortInt)                    :: idx

    if (p % fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = p % w

      ! Fast or thermal?
      if (p % E <= self % thermalEnergy) then
        idx = THERMAL_LEAK
      else
        idx = FAST_LEAK
      end if

      ! Score analog leakage
      call mem % score(histWgt, self % getMemAddress() + idx)

    end if

  end subroutine reportHist

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine closeCycle(self, end, mem)
    class(sixFactorClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    type(scoreMemory), intent(inout)     :: mem
    integer(longInt)                     :: addr
    real(defReal)                        :: fastFissProd, fastLeak, fastAbs,&
                                            thermalLeak, thermalAbs, thermalFuelAbs,&
                                            thermalFissProd, totalLoss, thermalLoss
    real(defReal)                        :: epsilon, fastPNL, p, thermalPNL, f, eta, keff

    if (mem % lastCycle()) then

      addr = self % getMemAddress()

      fastLeak = mem % getScore(addr + FAST_LEAK)
      fastAbs = mem % getScore(addr + FAST_ABS)
      thermalLeak = mem % getScore(addr + THERMAL_LEAK)
      thermalAbs = mem % getScore(addr + THERMAL_ABS)
      fastFissProd = mem % getScore(addr + FAST_PROD)
      thermalFissProd = mem % getScore(addr + THERMAL_PROD)
      thermalFuelAbs = mem % getScore(addr + THERMAL_FUEL_ABS)
      
      totalLoss = fastAbs + thermalAbs + fastLeak + thermalLeak
      thermalLoss = thermalAbs + thermalLeak

      epsilon = ONE + fastFissProd / thermalFissProd
      fastPNL = ONE - fastLeak / totalLoss
      p = thermalLoss / (fastAbs + thermalLoss)
      thermalPNL = thermalAbs / thermalLoss
      f = thermalFuelAbs / thermalAbs
      eta = thermalFissProd / thermalFuelAbs
      keff = epsilon * fastPNL * p * thermalPNL * f * eta

      ! Catch possible divisions by zero
      if (epsilon /= epsilon) epsilon = ZERO
      if (fastPNL /= fastPNL) fastPNL = ZERO
      if (p /= p) p = ZERO
      if (thermalPNL /= thermalPNL) thermalPNL = ZERO
      if (f /= f) f = ZERO
      if (eta /= eta) eta = ZERO
      if (keff /= keff) keff = ZERO

      call mem % accumulate(epsilon, addr + FAST_FISSION_FACTOR)
      call mem % accumulate(fastPNL, addr + FAST_NON_LEAKAGE)
      call mem % accumulate(p, addr + RESONANCE_ESCAPE_PROB)
      call mem % accumulate(thermalPNL, addr + THERMAL_NON_LEAKAGE)
      call mem % accumulate(f, addr + THERMAL_UTILISATION)
      call mem % accumulate(eta, addr + ETA_THERMAL)
      call mem % accumulate(keff, addr + K_EFFECTIVE)

    end if

  end subroutine closeCycle

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(sixFactorClerk), intent(in) :: self
    type(scoreMemory), intent(in)     :: mem
    real(defReal)                     :: epsilon, STDe, fastPNL, STDfPNL, p, STDp, &
                                         thermalPNL, STDtPNL, f, STDf, eta, STDeta, &
                                         keff, STDk
    character(MAX_COL)                :: buffer

    ! Get current factor estimates
    call mem % getResult(epsilon, STDe, self % getMemAddress() + FAST_FISSION_FACTOR)
    call mem % getResult(fastPNL, STDfPNL, self % getMemAddress() + FAST_NON_LEAKAGE)
    call mem % getResult(p, STDp, self % getMemAddress() + RESONANCE_ESCAPE_PROB)
    call mem % getResult(thermalPNL, STDtPNL, self % getMemAddress() + THERMAL_NON_LEAKAGE)
    call mem % getResult(f, STDf, self % getMemAddress() + THERMAL_UTILISATION)
    call mem % getResult(eta, STDeta, self % getMemAddress() + ETA_THERMAL)
    call mem % getResult(keff, STDk, self % getMemAddress() + K_EFFECTIVE)

    ! Print to console
    write (buffer, '(A,F8.5,A,F8.5)') 'Fast fission factor: ', epsilon, ' +/- ', STDe
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'Fast non-leakage probability: ', fastPNL, ' +/- ', STDfPNL
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'Resonance escape probability: ', p, ' +/- ', STDp
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'Thermal non-leakage probability: ', thermalPNL, ' +/- ', STDtPNL
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'Thermal utilisation: ', f, ' +/- ', STDf
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'Thermal eta: ', eta, ' +/- ', STDeta
    call statusMsg(buffer)
    write (buffer, '(A,F8.5,A,F8.5)') 'k-effective: ', keff, ' +/- ', STDk
    call statusMsg(buffer)

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(sixFactorClerk), intent(in) :: self
    class(outputFile), intent(inout)  :: outFile
    type(scoreMemory), intent(in)     :: mem
    character(nameLen)                :: name
    real(defReal)                     :: val, STD
    integer(longInt)                  :: addr

    call outFile % startBlock( self % getName())
    addr = self % getMemAddress()

    name = 'FAST_FISSION_FACTOR'
    call mem % getResult(val, STD, addr + FAST_FISSION_FACTOR)
    call outFile % printResult(val, STD, name)

    name = 'FAST_NON_LEAKAGE_P'
    call mem % getResult(val, STD, addr + FAST_NON_LEAKAGE)
    call outFile % printResult(val, STD, name)

    name = 'RESONANCE_ESCAPE_P'
    call mem % getResult(val, STD, addr + RESONANCE_ESCAPE_PROB)
    call outFile % printResult(val, STD, name)

    name = 'THERMAL_NON_LEAKAGE_P'
    call mem % getResult(val, STD, addr + THERMAL_NON_LEAKAGE)
    call outFile % printResult(val, STD, name)

    name = 'THERMAL_UTILISATION'
    call mem % getResult(val, STD, addr + THERMAL_UTILISATION)
    call outFile % printResult(val, STD, name)

    name = 'K_EFF_SIX_FACTOR'
    call mem % getResult(val, STD, addr + K_EFFECTIVE)
    call outFile % printResult(val, STD, name)

    call outFile % endBlock()

  end subroutine print

end module sixFactorClerk_class
