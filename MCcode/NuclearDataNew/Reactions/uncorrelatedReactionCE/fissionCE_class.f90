module fissionCE_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard
  use reactionHandle_inter,         only : reactionHandle
  use releaseLawENDF_inter,         only : releaseLawENDF
  use energyLawENDF_inter,          only : energyLawENDF
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  ! Data Structures
  use endfTable_class,              only : endfTable

  ! Factiories
  use releaseLawENDFfactory_func,   only : new_totalNu, new_delayedNu
  use energyLawENDFfactory_func,    only : new_energyLawENDF


  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: fissionCE_TptrCast

  !!
  !! Helper type to group together information for a delayed group
  !!
  !! Public Members:
  !!   lambda -> Decay constant of the precursor [1/s]
  !!   prob   -> Table with energy-dependent probability of this precursor
  !!   eLaw   -> Outgoing energy distribution from the precursor
  !!
  type, private :: precursor
    real(defReal)                    :: lambda
    type(endfTable)                  :: prob
    class(energyLawENDF),allocatable :: eLaw
  end type


  !!
  !! Contains information about Fission reaction based on ACE format
  !!
  !! Private Members:
  !!   nuBarTotal   -> Total (prompt & delayed) release table for incident energy [MeV]
  !!   eLawPrompt   -> Energy Law for the prompt fission neutrons
  !!   nuBarDelayed -> Delayed release table for incindent energy [MeV]
  !!   delayed      -> Information about Delayed emission precursors
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> construct and instance from ACE Card Data Deck
  !!
  type, public, extends(uncorrelatedReactionCE) :: fissionCE
    private
    class(releaseLawENDF),allocatable        :: nuBarTotal
    class(energyLawENDF), allocatable        :: eLawPrompt
    class(releaseLawENDF),allocatable        :: nuBarDelayed
    type(precursor),dimension(:),allocatable :: delayed

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: inCMframe
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleOut
    procedure :: probOf

    ! Type specific procedures
    procedure :: buildFromACE
  end type fissionCE

contains

  !!
  !! Initialsie
  !!
  !! See reactionHandle for details
  !!
  !! Errors:
  !!   fatalError if MT /= N_FISSION
  !!
  subroutine init(self, data, MT)
    class(fissioNCE), intent(inout) :: self
    class(dataDeck), intent(inout)  :: data
    integer(shortInt), intent(in)   :: MT
    character(100),parameter :: Here ='init (fissionCE_class.f90)'

    if( MT /= N_FISSION) then
      call fatalError(Here,'fissionCE suports only MT=18. Was given: '//numToChar(MT))
    end if

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
      type is (aceCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Fission CE cannot be build from '//data % myType())
    end select


  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fissionCE), intent(inout) :: self
    integer(shortInt)               :: i

    if(allocated(self % nuBarTotal)) then
      call self % nuBarTotal % kill()
      deallocate(self % nuBarTotal)
    end if

    if(allocated(self % eLawPrompt)) then
      call self % eLawPrompt % kill()
      deallocate(self % eLawPrompt)
    end if

    if(allocated(self % nuBarDelayed)) then
      call self % nuBarDelayed % kill()
      deallocate(self % nuBarDelayed)
    end if

    if(allocated(self % delayed)) then
      do i=1,size(self % delayed)
        call self % delayed(i) % prob % kill()

        if(allocated( self % delayed(i) % eLaw)) then
          call self % delayed(i) % eLaw % kill()
          deallocate(self % delayed(i) % eLaw)

        end if
      end do
      deallocate(self % delayed)
    end if

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(fissionCE), intent(in) :: self
    logical(defBool)             :: isIt

     isIt = .false.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function release(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = self % nuBarTotal % releaseAt(E)

  end function release

  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function releasePrompt(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = self % nuBarTotal % releaseAt(E)

  end function releasePrompt

  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand, lambda)
    class(fissionCE), intent(in)         :: self
    real(defReal), intent(out)           :: mu
    real(defReal), intent(out)           :: phi
    real(defReal), intent(out)           :: E_out
    real(defReal), intent(in)            :: E_in
    class(RNG), intent(inout)            :: rand
    real(defReal), intent(out), optional :: lambda

    ! Sample mu
    mu = TWO * rand % get() - ONE

    ! Sample Phi
    phi = TWO_PI * rand % get()

    ! Sample E_out
    E_out = self % eLawPrompt % sample(E_in, rand)

    ! Only prompt particles. Set delay
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)                :: mu
    real(defReal), intent(in)                :: phi
    real(defReal), intent(in)                :: E_out
    real(defReal), intent(in)                :: E_in
    real(defReal)                            :: prob

    if(abs(mu) <= ONE .and. E_out > ZERO .and. phi <= TWO_PI .and. phi >= ZERO) then
      prob = self % eLawPrompt % probabilityOf(E_out, E_in) / (TWO * TWO_PI)
    else
      prob = ZERO
    end if

  end function probOf

  !!
  !! Build fissionCE from ACE dataCard
  !!
  !! Args:
  !!   ACE [inout] -> ACE Card with the data
  !!
  !! Errors:
  !!   fatalError if there is no NU data (JXS(2) == 0)
  !!
  subroutine buildFromACE(self, ACE)
    class(fissionCE), intent(inout)            :: self
    type(aceCard), intent(inout)               :: ACE
    logical(defBool)                           :: withDelayed, onlyOneNu
    integer(shortInt)                          :: i, NR, N
    integer(shortInt),dimension(:),allocatable :: nrDat
    character(100),parameter :: Here = 'buildFromACE (fissionCE_class.f90)'

    call self % kill()

    ! Identify available data
    onlyOneNu   = ACE % hasNuTotal() .and. (.not.ACE % hasNuPrompt())
    withDelayed = ACE % hasNuDelayed()

    ! Detect unexpected data
    if(withDelayed .and. onlyOneNu) then
      call fatalError(Here, 'Prompt/Total Nu is given with delayed data. Which one is which? '//trim(ACE % ZAID))

    else if ( .not.ACE % hasNuTotal() .and. withDelayed) then
      call fatalError(Here, 'Has delayed neutron data but does not have total NuBar. WTF? '//trim(ACE % ZAID))

    end if

    ! Read basic data
    call new_totalNU(self % nuBarTotal, ACE)
    call new_energyLawENDF(self % eLawPrompt, ACE, N_FISSION)

    ! Read Delayed Data
    if(withDelayed) then
      ! Read Table
      call new_delayedNu(self % nuBarDelayed, ACE)

      ! Allocate data
      allocate (self % delayed(ACE % precursorGroups()))

      ! Detect missing data
      if(size(self % delayed) == 0) then
        call fatalError(Here, 'Has delayed neutrons but not precursors. WTF? '//trim(ACE % ZAID))
      end if

      ! Read Precursor data
      call ACE % setToPrecursors()
      do i=1,size(self % delayed)
        ! Read delay constant
        self % delayed(i) % lambda = ACE % readReal()
        nr = ACE % readInt()

        if(nr < 0) call fatalError(Here, 'NR < 0. WTF?')

        if(nr == 0) then ! Single interpolation region lin-lin
          N = ACE % readInt()
          associate ( dat => ACE % readRealArray(2*N))
            call self % delayed(i) % prob % init(dat(1:N), dat(N+1:2*N))
          end associate

        else ! Multiple Intrpolation regions
          nrDat = ACE % readIntArray(2*NR)
          N = ACE % readInt()
          associate (dat => ACE % readRealArray(2*N))
            call self % delayed(i) % prob % init(dat(1:N), dat(N+1:2*N), nrDat(1:NR), nrDat(NR+1:2*NR))
          end associate

        end if
      end do

      ! Read Energy distributions
      do i=1,size(self % delayed)
        call new_energyLawENDF(self % delayed(i) % eLaw, ACE, i, delayed = .true.)
      end do
    end if

  end subroutine buildFromACE

  !!
  !! Cast reactionHandle pointer to fissionCE pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of fissionCE type
  !!   Target points to source if source is fissionCE type
  !!
  pure function fissionCE_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(fissionCE), pointer                   :: ptr

    select type(source)
      type is(fissionCE)
        ptr => source

      class default
        ptr => null()

    end select

  end function fissionCE_TptrCast

end module fissionCE_class
