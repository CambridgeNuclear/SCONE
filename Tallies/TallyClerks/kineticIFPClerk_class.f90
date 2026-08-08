module kineticIFPClerk_class

  use numPrecision
  use tallyCodes
  use universalVariables,    only : MAX_COL
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError, numToChar
  use display_func,          only : statusMsg
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  implicit none
  private

  integer(shortInt), parameter :: N_PREC = 8       ! Fixed number of precursor groups. Some data
                                                   ! libraries will have only 6 groups and so
                                                   ! will have tallies of 0 in groups 7 and 8.
  integer(shortInt), parameter :: MEM_SIZE = N_PREC + 2 ! Total memory size
  
  ! Locations of different bins wrt memory address of the clerk
  integer(longInt), parameter  :: REMOV_TIME = 0 ,&  ! Address of removal time
                                  BETA_EFF = 1       ! Address of beta effective

  !!
  !! Iterated fission probability-based tally for kinetic parameters.
  !! These are beta-effective and the effective removal time.
  !! Meaningful results require that the calculation is run with a non-zero
  !! IFP generation size. Otherwise outputs will be zero.
  !! Note that removal time can be converted to generation time by division
  !! by keff.
  !!
  !! Private Members:
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type kineticIFPClerk;
  !! }
  !!
  type, public,extends(tallyClerk) :: kineticIFPClerk
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type kineticIFPClerk


contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(kineticIFPClerk), intent(inout) :: self
    class(dictionary), intent(in)         :: dict
    character(nameLen), intent(in)        :: name
    character(100),parameter :: Here = 'init (kineticIFPClerk.f90)'

    ! Needs no settings, just load name
    call self % setName(name)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(kineticIFPClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(kineticIFPClerk),intent(in)          :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(kineticIFPClerk), intent(in) :: self
    integer(shortInt)                  :: S

    S = MEM_SIZE

  end function getSize

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(kineticIFPClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    type(scoreMemory), intent(inout)      :: mem
    real(defReal)                         :: tau, beta
    real(defReal), dimension(N_PREC)      :: betaG
    integer(shortInt)                     :: g

    ! For each particle in the dungeon, obtain its ancestry.
    ! The IFP-weighted estimator for removal time is:
    ! For each particle:
    !   sum over ancestors(weight * lifetime) / sum over ancestors(weight)
    !
    ! The IFP-weighted estimator for beta-efffective is:
    ! For each particle:
    !   sum over ancestors(weight * delta_(born from prec)) / sum over ancestors(weight)
    ! 
    ! The beta estimator can be decomposed into groups by setting the delta = 1
    ! only if the neutron is born from a particular delayed group
    !
    ! These quantities are computed in the dungeon, averaged over particles, 
    ! to prevent copying many arrays into the tally.

    call end % getIFPValues(N_PREC, tau, beta, betaG)
    call mem % accumulate(tau, self % getMemAddress() + REMOV_TIME)
    call mem % accumulate(beta, self % getMemAddress() + BETA_EFF)
    do g = 1, N_PREC
      call mem % accumulate(betaG(g), self % getMemAddress() + BETA_EFF + g)
    end do

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(kineticIFPClerk), intent(in)  :: self
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: tau, STDtau, beta, STDbeta
    character(MAX_COL)                  :: buffer

    call mem % getResult(tau, STDtau, self % getMemAddress() + REMOV_TIME)
    call mem % getResult(beta, STDbeta, self % getMemAddress() + BETA_EFF)

    ! Print estimates to a console
    write(buffer, '(A,ES12.5,A,ES12.5)') 'Removal time (IFP): ', tau, ' +/- ', STDtau
    call statusMsg(buffer)
    write(buffer, '(A,ES12.5,A,ES12.5)') 'Beta-effective (IFP): ', beta, ' +/- ', STDbeta
    call statusMsg(buffer)

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(kineticIFPClerk), intent(in) :: self
    class(outputFile), intent(inout)   :: outFile
    type(scoreMemory), intent(in)      :: mem
    real(defReal)                      :: tau, beta, STDtau, STDbeta
    real(defReal), dimension(N_PREC)   :: betaG, STDbetaG
    integer(shortInt)                  :: g
    character(nameLen)                 :: name

    ! Get result values
    call mem % getResult(tau, STDtau, self % getMemAddress() + REMOV_TIME)
    call mem % getResult(beta, STDbeta, self % getMemAddress() + BETA_EFF)

    do g = 1, N_PREC
      call mem % getResult(betaG(g), STDbetaG(g), self % getMemAddress() + BETA_EFF + g)
    end do

    ! Print to output file
    call outFile % startBlock(self % getName())
    name = 'Removal_time'
    call outFile % printResult(tau, STDtau, name)
    name = 'Beta_effective'
    call outFile % printResult(beta, STDbeta, name)
    do g = 1, N_PREC
      name = 'Beta_g'//numToChar(g)
      call outFile % printResult(betaG(g), STDbetaG(g), name)
    end do
    call outFile % endBlock()

  end subroutine print

end module kineticIFPClerk_class
