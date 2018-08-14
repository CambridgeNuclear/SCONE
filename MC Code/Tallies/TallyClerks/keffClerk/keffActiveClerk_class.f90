module keffActiveClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use genericProcedures,          only : fatalError, charCmp
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use keffClerk_inter,            only : keffClerk
  use tallyEstimator_class,       only : tallyScore, tallyCounter
  use transportNuclearData_inter, only : transportNuclearData
  use outputFile_class,           only : outputFile

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'keffActiveClerk'

  interface keffActiveClerk
    module procedure new_keffActiveClerk
  end interface

  type, public,extends(keffClerk) :: keffActiveClerk
    private
    character(nameLen)   :: name
    integer(shortInt)    :: cycleCount = 0

    type(tallyScore)     :: impProd     ! Implicit neutron production
    type(tallyScore)     :: scatterProd ! Analog estimate of scattering production by (N,XN)
    type(tallyScore)     :: impAbs      ! Implicit neutron absorbtion
    type(tallyScore)     :: anaLeak     ! Analog neutron leakage

    type(tallyCounter)   :: k_analog
    type(tallyCounter)   :: k_imp

    real(defReal)        :: startWgt
    real(defReal)        :: targetSD = 0.0

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    procedure :: keff
    procedure :: print

    ! Overwrite report procedures
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
  end type keffActiveClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(keffActiveClerk), intent(in)         :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, outColl_CODE, cycleStart_CODE, cycleEnd_CODE, hist_CODE]

  end function validReports

  !!
  !! Display progress current estimate of k-eff with STD
  !!
  subroutine display(self)
    class(keffActiveClerk), intent(in) :: self
    real(defReal)                      :: k_imp, k_analog, STD_imp, STD_analog

    ! Obtain current estimates of k analog and implicit
    call self % k_imp % getEstimate(k_imp, STD_imp, self % cycleCount)
    call self % k_analog % getEstimate(k_analog, STD_analog, self % cycleCount)

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'k-eff (implicit): ', k_imp, ' +/- ', STD_imp
    print '(A,F8.5,A,F8.5)', 'k-eff (analog): ',  k_analog, ' +/- ', STD_analog

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(keffActiveClerk), intent(in) :: self
    logical(defBool)                   :: isIt
    real(defReal)                      :: k, SD

    call self % k_imp % getEstimate(k,SD,self % cycleCount)

    isIt = (SD < self % targetSD)

  end function isConverged


  !!
  !! Process collision report
  !!
  subroutine reportInColl(self,p)
    class(keffActiveClerk), intent(inout) :: self
    class(particle), intent(in)           :: p
    type(xsMacroSet_ptr)                  :: XSs
    real(defReal)                         :: totalXS, nuFissXS, absXS, flux
    real(defReal)                         :: s1, s2
    character(100), parameter  :: Here = 'reportInColl (keffActiveClerk_class.f90)'

    ! Obtain XSs
    ! Check if it dynamic type is supported
    ! If it is obtain macroscopic XSs
    ! It it isn't throw error
    associate (xsData => p % xsData)
      select type(xsData)
        class is (transportNuclearData)
          call xsData % getMatMacroXS(XSs, p, p % matIdx())

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    totalXS  = XSs % totalXS()
    nuFissXS = XSs % nuFissionXS()
    absXS    = XSs % captureXS() + XSs % fissionXS()

    ! Calculate flux and scores
    flux = p % w / totalXS

    s1 = nuFissXS * flux
    s2 = absXS * flux

    ! Add scores to counters
    call self % impProd % add(s1)
    call self % impAbs  % add(s2)

  end subroutine reportInColl

  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self,pre,post,MT,muL)
    class(keffActiveClerk), intent(inout) :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
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
    if (mult /= ZERO) call self % scatterProd % add(pre % wgt * mult)

  end subroutine reportOutColl

  !!
  !! Process history report
  !!
  subroutine reportHist(self,pre,post,fate)
    class(keffActiveClerk), intent(inout):: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate
    real(defReal)   :: histWgt


    if( fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = pre % wgt

      ! Score analog leakage
      call self % anaLeak % add(histWgt)

    end if

  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(keffActiveClerk), intent(inout):: self
    class(particleDungeon), intent(in)   :: start

    self % startWgt  = start % popWeight()

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(keffActiveClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    real(defReal)                         :: endWgt, k_est
    real(defReal)                         :: nuFiss, absorb, leakage, scatterMul, k_cycle

    ! Obtain end of cycle weight and k value used to change fission site generation rate
    endWgt = end % popWeight()
    k_cycle = end % k_eff

    ! Calculate and score analog estimate of k-eff
    k_est =  endWgt / self % startWgt * k_cycle
    call self % k_analog % addEstimate(k_est)

    ! Calculate and score implicit estimate of k_eff
    nuFiss     = self % impProd % get()
    absorb     = self % impAbs % get()
    leakage    = self % anaLeak % get()
    scatterMul = self % scatterProd % get()

    k_est = nuFiss / (absorb + leakage - scatterMul )

    call self % k_imp % addEstimate(k_est)

    ! Reset score counters
    call self % impProd % reset()
    call self % impAbs % reset()
    call self % anaLeak % reset()
    call self % scatterProd % reset()

    ! Increas counter of cycles
    self % cycleCount = self % cycleCount + 1

  end subroutine reportCycleEnd

  !!
  !! Write contents of the keffActive Clerk to output file
  !!
  subroutine print(self,outFile)
    class(keffActiveClerk), intent(in) :: self
    class(outputFile), intent(inout)   :: outFile
    real(defReal)                      :: val, std
    character(nameLen)                 :: name
    call outFile % startBlock(self % name)


    !call self % k_analog % getEstimate(k_analog, STD_analog, self % cycleCount)

    ! Print k implicit
    call self % k_imp % getEstimate(val, std, self % cycleCount)
    name = 'k_Imp'
    call outFile % printResult(val,std,name)

    ! Print k implicit
    call self % k_analog % getEstimate(val, std, self % cycleCount)
    name = 'k_Analog'
    call outFile % printResult(val,std,name)

    call outFile % endBlock()


  end subroutine print

  !!
  !! Initialise keffActiveClerk from dictionary
  !!
  subroutine init(self,dict,name)
    class(keffActiveClerk),intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(nameLen)                   :: type
    character(100),parameter :: Here ='init (keffActiveClerk_class.f90)'

    ! Assign name
    self % name = name

    call dict % get(type,'type')

    ! Check that class description matches class name
!    if ( .not.charCmp(CLASS_NAME,type) ) then
!      call fatalError(Here, 'Type : ' // type // ' is different form class name ' // CLASS_NAME )
!
!    end if

    call dict % getOrDefault(type,'trigger','no')

    ! Read convergance target
    if( charCmp(type,'yes')) then
      call dict % get(self % targetSD,'SDtarget')

    end if

  end subroutine init


  !!
  !! Return estimate of k-eff
  !!
  function keff(self) result(k)
    class(keffActiveClerk), intent(in) :: self
    real(defReal)                      :: k

    call self % k_imp % getEstimate(k, self % cycleCount)

  end function keff

  !!
  !! keffActiveClerk constructor function
  !!
  function new_keffActiveClerk(dict, name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(keffActiveClerk)          :: new

    call new % init(dict, name)

  end function new_keffActiveClerk

end module keffActiveClerk_class
