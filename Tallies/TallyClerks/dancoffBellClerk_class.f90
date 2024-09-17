module dancoffBellClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables,         only : MAX_COL
  use genericProcedures,          only : fatalError, hasDuplicates
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile
  use intMap_class,               only : intMap

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill
  use tallyResult_class,          only : tallyResult
  use energyFilter_class,         only : energyFilter

  ! Nuclear Data
  use materialMenu_mod,           only : mm_matIdx => matIdx
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !! Local parameters. Note that OUSIDE is not arbitraty. Choosen to match invalid idx
  !! given by a tallyMap
  integer(shortInt), parameter :: FUEL      = -2, &
                                  MODERATOR = -3, &
                                  OUTSIDE   = 0

  !! Offsets for diffrent bins
  integer(longInt), parameter :: ESC_PROB_TOTXS = 0, &
                                 STAY_PROB      = 1, &
                                 D_EFF          = 2


  !!
  !! Tally clerk to calculate a single-term rational approximation for a lattice
  !!
  !! Calculates a term in a single-term rational approximation by totalXS moments
  !! of escape probabbility. Calculates: D_eff = P_e * SIGMA_t / (1-P_e)
  !!
  !! In MC calculatiuon D_Effcan be calculated per batch as:
  !! D_eff = wgt_e * SIGMA_t / wgt_fuel
  !!
  !! Where wgt_e is weight that escapes fuel regions and wgt_mod is weight that remains in fuel.
  !!
  !! Private members:
  !!   materialSet -> IntMap that stores relevant matIdxs and designations (FUEL or MODERATOR)
  !!   filter      -> Energy filter
  !!
  !! Interface:
  !!   tallyClerk interface
  !!
  !! Sample Input Dictionary:
  !!
  !! myClerk {
  !!   type dancoffBellClerk;
  !!   # Etop  10;   #           // Top energy boundary [MeV] (def: 20)
  !!   # Elow  0.06; #           // Bottom energy boundary [MeV] (def: 0)
  !!   fuelMat (m1 m2 m3);       // List of fuel material names
  !!   modMat  (m4 m5 m98);      // List of moderator material names
  !!  }
  !!
  type, public,extends(tallyClerk) :: dancoffBellClerk
    private
    type(intMap)         :: materialSet
    type(energyFilter)   :: filter
  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportTrans
    procedure  :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill
  end type dancoffBellClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(dancoffBellClerk), intent(inout)      :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    real(defReal)                               :: Emax, Emin
    character(nameLen),dimension(:),allocatable :: fuelNames, modNames
    integer(shortInt),dimension(:), allocatable :: fuelIdx, modIdx
    integer(shortInt)                           :: i
    character(100), parameter :: Here ='init (dancoffBellClerk_class.f90)'

    ! Load name
    call self % setName(name)

    ! Initialise filter
    call dict % getOrDefault(Emin,'Elow', ZERO)
    call dict % getOrDefault(Emax,'Etop', 20.0_defReal)
    call self % filter % build(Emin, Emax)

    ! Get fuel and moderator materials
    call dict % get(fuelNames, 'fuelMat')
    call dict % get(modNames, 'modMat')

    ! Convert to indexes
    allocate(fuelIdx(size(fuelNames)))
    allocate(modIdx(size(modNames)))

    do i = 1,size(fuelNames)
      fuelIdx(i) = mm_matIdx(fuelNames(i))
    end do

    do i =1,size(modNames)
      modIdx(i) = mm_matIdx(modNames(i))
    end do

    ! Check for overlaps
    ! Check for overlap
    if (hasDuplicates([modIdx, fuelIdx])) then
      call fatalError(Here, 'Same materials are defined as fuel and moderator')
    end if

    ! Create map of matIdx -> matType
    call self % materialSet % init (size(modIdx) + size(fuelIdx))

    ! Load Fuel
    do i = 1,size(fuelIdx)
      call self % materialSet % add(fuelIdx(i), FUEL)
    end do

    ! Load Moderator
    do i = 1,size(modIdx)
      call self % materialSet % add(modIdx(i), MODERATOR)
    end do

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(dancoffBellClerk),intent(in)         :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [trans_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(dancoffBellClerk), intent(in) :: self
    integer(shortInt)                   :: S

    S = 3

  end function getSize

  !!
  !! Process transition report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTrans(self, p, xsData, mem)
    class(dancoffBellClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    class(nuclearDatabase),intent(inout)   :: xsData
    type(scoreMemory), intent(inout)       :: mem
    real(defReal)                          :: SigmaTot
    integer(shortInt)                      :: T_end, T_start
    real(defReal)                          :: w_end
    type(particleState)                    :: state

    character(100),parameter :: Here = 'reportTrans (dancoffBellClerk_class.f90)'

    ! Find start material type; Exit if not fuel
    T_start = self % materialSet % getOrDefault(p % preTransition % matIdx, OUTSIDE)
    if (T_start /= FUEL) return

    ! Exit if outside energy range
    state = p
    if (.not.self % filter % isPass(state)) return

    ! Find end material type; Exit if not fuel or moderator
    T_end = self % materialSet % getOrDefault(p % matIdx(), OUTSIDE)
    if (T_end == OUTSIDE) return

    ! Obtain starting and ending weights
    w_end   = p % w

    ! Add to approperiate bins
    select case (T_end)
      case (MODERATOR)
        ! Get XS
        SigmaTot = xsData % getTotalMatXS(p, p % preTransition % matIdx)

        call mem % score(w_end * SigmaTot, self % getMemAddress() + ESC_PROB_TOTXS)

      case (FUEL)
        call mem % score(w_end, self % getMemAddress() + STAY_PROB)

      case default
        call fatalError(Here, 'WTF? Impossible state')

    end select

  end subroutine reportTrans

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(dancoffBellClerk), intent(inout)  :: self
    class(particleDungeon), intent(in)      :: end
    type(scoreMemory), intent(inout)        :: mem
    real(defReal)                           :: escSigmaT, fuelWgt

    if (mem % lastCycle()) then
      escSigmaT = mem % getScore(self % getMemAddress() + ESC_PROB_TOTXS)
      fuelWgt   = mem % getScore(self % getMemAddress() + STAY_PROB)
      call mem % accumulate(escSigmaT / fuelWgt, self % getMemAddress() + D_EFF)
    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(dancoffBellClerk), intent(in) :: self
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: mean, STD
    character(MAX_COL)                  :: buffer

    call mem % getResult(mean, STD, self % getMemAddress() + D_EFF)

    ! Print to console
    write (buffer, '(A,ES15.5,A,ES15.5)') 'Dancoff-Bell: ', mean, ' +/- ', STD
    call statusMsg(buffer)

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(dancoffBellClerk), intent(in) :: self
    class(outputFile), intent(inout)    :: outFile
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: mean, STD
    character(nameLen)                  :: name

    ! Get result
    call mem % getResult(mean, STD, self % getMemAddress() + D_EFF)

    ! Print result
    call outFile % startBlock( self % getName())

    name = 'DanBellFactor'
    call outFile % printResult(mean, STD, name)

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  !! See tallyClerk_inter for details
  !!
  elemental subroutine kill(self)
    class(dancoffBellClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    call self % materialSet % kill()
    ! Kill filter whan available

  end subroutine kill


end module dancoffBellClerk_class
