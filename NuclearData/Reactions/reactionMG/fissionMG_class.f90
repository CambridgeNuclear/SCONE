module fissionMG_class

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, numToChar
  use RNG_class,             only : RNG
  use reactionHandle_inter,  only : reactionHandle
  use reactionMG_inter,      only : reactionMG
  use dataDeck_inter,        only : dataDeck
  use dictDeck_class,        only : dictDeck
  use dictionary_class,      only : dictionary

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: fissionMG_TptrCast

  !!
  !! A special type of MG reaction that contains data related to fission
  !!
  !! At the moment it has no information about delayed emissions, shall be introduced
  !! at a latter date.
  !! Fission spectrum is independent of the energy group
  !!
  !! Public members:
  !!   data  -> data for fissionMG [energyGroup, reaction (nu or chi)]
  !!   kappa -> Energy per fission. 200 MeV by default.
  !!
  !! Interface:
  !!   reactionMG interface
  !!   buildFromDict -> builds fissionMG from a SCONE dictionary
  !!   getKappa      -> obtains the energy per fission in MeV
  !!
  type, public, extends(reactionMG) :: fissionMG
    real(defReal),dimension(:,:),allocatable :: data
    real(defReal) :: kappa = 200.0_defReal
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleDelayRate
    procedure :: sampleOut

    ! Local procedures
    procedure :: buildFromDict
    procedure :: getKappa

  end type fissionMG

  !!
  !! Reaction indices
  !!
  integer(shortInt),parameter :: NU_DAT = 1, CHI_DAT = 2

contains

  !!
  !! Initialise fissionMG from data deck
  !!
  !! See reactionHandle for details.
  !!
  !! Errors:
  !!   Returns a fatalError if any data deck diffrent from SCONE dictionary is given
  !!   Returns fatalError if MT is diffrent from macroFission
  !!
  subroutine init(self, data, MT)
    class(fissionMG), intent(inout) :: self
    class(dataDeck), intent(inout)  :: data
    integer(shortInt), intent(in)   :: MT
    character(100), parameter :: Here = 'init (fissionMG_class.f90)'

    ! Verify that MT is OK
    if (MT /= macroFission) then
      call fatalError(Here,'Unsupported MT number: '//numToChar(MT)//' only macroFission is allowed')
    end if

    ! Select approperiate build procedure for a
    select type(data)
      type is(dictDeck)
        call self % buildFromDict(data % dict)

      class default
        call fatalError(Here, 'fissionMG cannot be build from: '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fissionMG), intent(inout) :: self

    if(allocated(self % data)) deallocate(self % data)

  end subroutine kill

  !!
  !! Returns number of particles produced on average by the fission
  !!
  !! See reactionMG documentation for details
  !!
  pure function release(self, G) result(N)
    class(fissionMG), intent(in)  :: self
    integer(shortInt), intent(in) :: G
    real(defReal)                 :: N

    if( G > 0 .and. G <= size(self % data, 1)) then
      N = self % data(G, NU_DAT)

    else
      N = ZERO

    end if

  end function release

  !!
  !! Returns number of particles produced instantly on average by the fission
  !!
  !! See reactionMG documentation for details
  !!
  pure function releasePrompt(self, G) result(N)
    class(fissionMG), intent(in)  :: self
    integer(shortInt), intent(in) :: G
    real(defReal)                 :: N

    N = self % release(G)

  end function releasePrompt

  !!
  !! Returns number of particles produced with delay on average by the fission
  !!
  !! See reactionMG documentation for details
  !!
  !! NOTE:
  !!   For now returns always 0.0 !
  !!
  pure function releaseDelayed(self, G) result(N)
    class(fissionMG), intent(in)  :: self
    integer(shortInt), intent(in) :: G
    real(defReal)                 :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample the delay rate for the delayed particle
  !!
  !! See reactionMG documentation for details
  !!
  !! NOTE:
  !!   For now always returns 0.0 !
  !!
  function sampleDelayRate(self, G, rand) result(lambda)
    class(fissionMG), intent(in) :: self
    integer(shortInt), intent(in):: G
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: lambda

    lambda = ZERO

  end function sampleDelayRate

  !!
  !! Sample outgoing particle
  !!
  !! See reactionMG documentation for details
  !!
  subroutine sampleOut(self, mu, phi, G_out, G_in, rand)
    class(fissionMG), intent(in)   :: self
    real(defReal), intent(out)     :: mu
    real(defReal), intent(out)     :: phi
    integer(shortInt), intent(out) :: G_out
    integer(shortInt), intent(in)  :: G_in
    class(RNG), intent(inout)      :: rand
    real(defReal)                  :: rem
    character(100),parameter :: Here = 'sampleOut (fissionMG_class.f90)'

    ! Sample mu and phi -> isotropic
    mu = TWO * rand % get() - ONE
    phi = TWO_PI * rand % get()

    ! Sample G_out
    rem = rand % get()
    do G_out=1,size(self % data(:,CHI_DAT))
      rem = rem - self % data(G_out, CHI_DAT)
      if(rem < ZERO) return
    end do

    call fatalError(Here,'WTF? Sampling failed. Unnormalised CHI or rand above 1?!')

  end subroutine sampleOut

  !!
  !! Builds fissionMG from SCONE dictionary
  !!
  !! Args:
  !!   dict [in] -> dictionary that contains data
  !!
  !! Errors:
  !!   FatalError if required data is not present in the dictionary
  !!   FatalError if chi is not normalised to 1.0 within FP_REL_TOL
  !!
  subroutine buildFromDict(self, dict)
    class(fissionMG), intent(inout)        :: self
    class(dictionary), intent(in)          :: dict
    integer(shortInt)                      :: nG
    real(defReal)                          :: S
    real(defReal),dimension(:),allocatable :: temp
    character(100),parameter :: Here = 'buildFromDict (fissionMG_class.f90)'

    ! Get number of groups
    call dict % get(nG, 'numberOfGroups')
    
    ! Get energy per fission
    call dict % getOrDefault(self % kappa, 'kappa', 200.0_defReal)

    ! Allocate space
    allocate(self % data(nG, 2))

    ! Get nu
    call dict % get(temp, 'nu')
    if(size(temp) /= ng) then
      call fatalError(Here, 'Invalid number of values of nu. Given: '// numToChar(size(temp)) // &
                            ' Expected: ' // numToChar(nG))
    end if
    self % data(:,NU_DAT) = temp

    ! Get Chi
    call dict % get(temp, 'chi')
    if(size(temp) /= ng) then
      call fatalError(Here, 'Invalid number of values of chi. Given: '// numToChar(size(temp)) // &
                            ' Expected: ' // numToChar(nG))
    end if
    self % data(:,CHI_DAT) = temp

    ! Check normalisation of chi
    S = sum(self % data(:,CHI_DAT))
    if( abs(S-ONE) > 0.01 * FP_REL_TOL) then
      print *,'Chi is not normalised. Relative error wrt ONE:'// numToChar(abs(S-ONE))//&
              ' The normalisation has been adjusted automatically'
      self % data(:,CHI_DAT) = self % data(:,CHI_DAT) / S
    end if

  end subroutine buildFromDict

  !!
  !! Get the energy release from fission
  !!
  pure function getKappa(self) result(kappa)
    class(fissionMG), intent(in) :: self
    real(defReal)                :: kappa

    kappa = self % kappa

  end function getKappa


  !!
  !! Cast reactionHandle pointer to fissionMG pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of fissionMG type
  !!   Target points to source if source is fissionMG type
  !!
  pure function fissionMG_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(fissionMG), pointer                   :: ptr

    select type(source)
      type is(fissionMG)
        ptr => source

      class default
        ptr => null()
    end select

  end function fissionMG_TptrCast


end module fissionMG_class
